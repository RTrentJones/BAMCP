"""E2E tests for viewer chrome: legend, no-reference banner, hidden-reads indicator,
UI state sync, base-zoom strand/MAPQ cues and loadData defensiveness."""

import pytest
from playwright.sync_api import Page

from bamcp.resources import get_viewer_html
from tests.e2e.conftest import sample_pixel_at_genomic

pytestmark = pytest.mark.e2e


@pytest.fixture
def viewer_page(page: Page):
    page.set_viewport_size({"width": 1024, "height": 800})
    page.set_content(get_viewer_html())
    page.wait_for_load_state("domcontentloaded")
    page.wait_for_function("() => typeof window.viewer !== 'undefined'", timeout=5000)
    return page


# Mirrors bamcp.eval.renderer._HYDRATE_JS: state is mutated directly, then resize().
HYDRATE_JS = """
(({payload, mode}) => {
    viewer.state.loadData(payload);
    viewer.state.settings.displayMode = mode;
    viewer.state.resortAndRepack();
    viewer.renderVariantTable();
    viewer.renderer.resize();
})
"""


def _read(name, pos, length=50, mapq=60, reverse=False, seq=None):
    return {
        "name": name,
        "sequence": seq if seq is not None else ("ACGTACGTAC" * 20)[:length],
        "cigar": f"{length}M",
        "position": pos,
        "end_position": pos + length,
        "mapping_quality": mapq,
        "is_reverse": reverse,
        "mismatches": [],
    }


def _payload(reads, start=100, end=160, **extra):
    span = end - start
    return {
        "contig": "chr1",
        "start": start,
        "end": end,
        "reads": reads,
        "coverage": [len(reads)] * span,
        "variants": [],
        "reference_sequence": ("ACGTACGTAC" * (span // 10 + 1))[:span],
        **extra,
    }


def pixel_in_cell(page: Page, genomic_pos: int, y: int, dx: int = 2) -> dict:
    """Sample a pixel `dx` px inside the base cell at genomic_pos (avoids cell edges)."""
    return page.evaluate(
        """([pos, y, dx]) => {
            const c = document.getElementById('reads-canvas');
            const scale = c.width / (viewer.state.viewport.end - viewer.state.viewport.start);
            const x = Math.floor((pos - viewer.state.viewport.start) * scale) + dx;
            const p = c.getContext('2d').getImageData(x, y, 1, 1).data;
            return {r: p[0], g: p[1], b: p[2], a: p[3]};
        }""",
        [genomic_pos, y, dx],
    )


def hydrate(page: Page, payload: dict, mode: str = "compact") -> None:
    page.evaluate(HYDRATE_JS, {"payload": payload, "mode": mode})
    page.wait_for_timeout(100)


def _lum(px: dict) -> float:
    """WCAG relative luminance of a sampled pixel."""

    def chan(v: int) -> float:
        c = v / 255
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4

    return 0.2126 * chan(px["r"]) + 0.7152 * chan(px["g"]) + 0.0722 * chan(px["b"])


def _contrast(a: dict, b: dict) -> float:
    la, lb = sorted((_lum(a), _lum(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


class TestLegend:
    def test_legend_outside_viewer_and_explains_encodings(self, viewer_page: Page):
        inside = viewer_page.evaluate("() => !!document.querySelector('#viewer #legend')")
        assert not inside, "legend must stay outside #viewer (eval clips to #viewer)"
        text = viewer_page.locator("#legend").inner_text()
        for label in (
            "Bases",
            "Strand",
            "Forward",
            "Reverse",
            "MAPQ",
            "Mismatch",
            "Insertion",
            "Deletion",
            "Coverage",
        ):
            assert label in text, f"legend missing {label!r}"

    def test_dv_legend_only_in_dv_mode(self, viewer_page: Page):
        hydrate(viewer_page, _payload([_read("r", 100)]), "compact")
        assert not viewer_page.locator(".legend-dv-composite").is_visible()
        hydrate(viewer_page, _payload([_read("r", 100)]), "dv-composite")
        assert viewer_page.locator(".legend-dv-composite").is_visible()


class TestReferenceBanner:
    def test_banner_shown_without_reference(self, viewer_page: Page):
        payload = _payload(
            [_read("r", 100)],
            applied_filters={"snv_detection": "disabled: no reference", "indel_detection": False},
        )
        hydrate(viewer_page, payload)
        banner = viewer_page.locator("#reference-banner")
        assert banner.is_visible()
        text = banner.inner_text()
        assert "No reference genome" in text
        assert "demo:HG002" in text

    def test_banner_hidden_when_detection_enabled(self, viewer_page: Page):
        hydrate(
            viewer_page, _payload([_read("r", 100)], applied_filters={"snv_detection": "enabled"})
        )
        assert not viewer_page.locator("#reference-banner").is_visible()

    def test_banner_hidden_when_fields_absent(self, viewer_page: Page):
        hydrate(viewer_page, _payload([_read("r", 100)]))
        assert not viewer_page.locator("#reference-banner").is_visible()

    def test_banner_clears_when_reference_arrives(self, viewer_page: Page):
        hydrate(
            viewer_page,
            _payload(
                [_read("r", 100)], applied_filters={"snv_detection": "disabled: no reference"}
            ),
        )
        assert viewer_page.locator("#reference-banner").is_visible()
        hydrate(
            viewer_page,
            _payload(
                [_read("r", 100)],
                applied_filters={"snv_detection": "enabled", "reference_source": "auto"},
            ),
        )
        assert not viewer_page.locator("#reference-banner").is_visible()


class TestHiddenReadsIndicator:
    def _deep_payload(self, n=120):
        # Every read overlaps every other → one read per packed row, n rows.
        return _payload([_read(f"r{i}", 100 + (i % 5)) for i in range(n)], start=0, end=400)

    def test_indicator_counts_reads_below_the_fold(self, viewer_page: Page):
        hydrate(viewer_page, self._deep_payload(), "expanded")
        ind = viewer_page.locator("#hidden-reads-indicator")
        assert ind.is_visible()
        count = int(ind.get_attribute("data-count"))
        rows = viewer_page.evaluate("() => viewer.state.packedRows.length")
        assert 0 < count < rows
        assert "reads not shown" in ind.inner_text()
        assert str(count) in ind.inner_text()

    def test_indicator_is_blind_hideable_inside_viewer(self, viewer_page: Page):
        hydrate(viewer_page, self._deep_payload(), "expanded")
        meta = viewer_page.evaluate(
            """() => {
                const el = document.getElementById('hidden-reads-indicator');
                return {
                    inViewer: !!el.closest('#viewer'),
                    blind: el.hasAttribute('data-blind-hide'),
                };
            }"""
        )
        assert meta == {"inViewer": True, "blind": True}

    def test_indicator_pinned_to_visible_bottom(self, viewer_page: Page):
        hydrate(viewer_page, self._deep_payload(), "expanded")
        box = viewer_page.evaluate(
            """() => {
                const v = document.getElementById('viewer').getBoundingClientRect();
                const i = document.getElementById('hidden-reads-indicator').getBoundingClientRect();
                return {vBottom: v.bottom, iBottom: i.bottom};
            }"""
        )
        assert abs(box["vBottom"] - box["iBottom"]) <= 2

    def test_indicator_hides_after_scrolling_to_bottom(self, viewer_page: Page):
        hydrate(viewer_page, self._deep_payload(), "expanded")
        viewer_page.evaluate(
            "() => { const v = document.getElementById('viewer'); v.scrollTop = v.scrollHeight; }"
        )
        viewer_page.wait_for_timeout(100)
        assert not viewer_page.locator("#hidden-reads-indicator").is_visible()

    def test_indicator_hidden_for_shallow_pileup(self, viewer_page: Page):
        hydrate(viewer_page, _payload([_read("a", 100), _read("b", 105, reverse=True)]))
        assert not viewer_page.locator("#hidden-reads-indicator").is_visible()


class TestUiStateSync:
    def test_mode_select_reflects_render_mode(self, viewer_page: Page):
        hydrate(viewer_page, _payload([_read("r", 100)]), "expanded")
        assert viewer_page.locator("#display-mode-select").input_value() == "expanded"
        hydrate(viewer_page, _payload([_read("r", 100)]), "dv-composite")
        assert viewer_page.locator("#display-mode-select").input_value() == "dv-composite"

    def test_region_box_shows_one_based_region_after_hydrate(self, viewer_page: Page):
        hydrate(viewer_page, _payload([_read("r", 100)], start=1039, end=1100))
        # Internal 0-based half-open [1039, 1100) is samtools chr1:1040-1100.
        assert viewer_page.locator("#region-input").input_value() == "chr1:1040-1100"

    def test_region_box_follows_zoom(self, viewer_page: Page):
        hydrate(viewer_page, _payload([_read("r", 100)], start=0, end=400))
        viewer_page.click("#zoom-in")
        viewer_page.wait_for_timeout(50)
        vp = viewer_page.evaluate("() => viewer.state.viewport")
        expected = f"chr1:{vp['start'] + 1}-{vp['end']}"
        assert viewer_page.locator("#region-input").input_value() == expected


class TestLoadDataDefensiveness:
    def test_missing_variants_does_not_throw(self, viewer_page: Page):
        errors = []
        viewer_page.on("pageerror", lambda e: errors.append(str(e)))
        payload = _payload([_read("r", 100)])
        del payload["variants"]
        hydrate(viewer_page, payload)
        assert errors == []
        assert viewer_page.evaluate("() => viewer.state.data.variants.length") == 0
        assert viewer_page.locator("#variant-table tr").count() == 0

    def test_missing_reads_and_coverage_does_not_throw(self, viewer_page: Page):
        errors = []
        viewer_page.on("pageerror", lambda e: errors.append(str(e)))
        hydrate(viewer_page, {"contig": "chr1", "start": 0, "end": 100})
        assert errors == []
        assert viewer_page.evaluate("() => viewer.state.packedRows.length") == 0


class TestBaseZoomCues:
    """At base zoom bases cover the read body; strand + MAPQ must stay visible."""

    # compact: 12px rows, 2px gap; strand stripe occupies the bottom 2px (y=10..11).
    STRIPE_Y = 11
    BASE_Y = 4

    def test_strand_stripe_colours_at_base_zoom(self, viewer_page: Page):
        payload = _payload([_read("fwd", 100), _read("rev", 105, reverse=True)])
        hydrate(viewer_page, payload)
        assert (
            viewer_page.evaluate(
                "() => document.getElementById('reads-canvas').width / "
                "(viewer.state.viewport.end - viewer.state.viewport.start)"
            )
            >= 10
        )
        fwd = sample_pixel_at_genomic(viewer_page, "reads-canvas", 120, self.STRIPE_Y)
        rev = sample_pixel_at_genomic(viewer_page, "reads-canvas", 120, 14 + self.STRIPE_Y)
        assert fwd["b"] > fwd["r"] + 60, f"forward stripe should be blue: {fwd}"
        assert rev["r"] > rev["b"] + 60, f"reverse stripe should be red: {rev}"

    def test_mapq_fade_at_base_zoom(self, viewer_page: Page):
        # The fade is an opaque blend toward white (theme-independent contrast), so a
        # low-MAPQ read's matching-base tiles are much lighter, not transparent.
        payload = _payload([_read("hi", 100, mapq=60), _read("lo", 105, mapq=6)])
        hydrate(viewer_page, payload)
        hi = pixel_in_cell(viewer_page, 121, self.BASE_Y)
        lo = pixel_in_cell(viewer_page, 121, 14 + self.BASE_Y)
        assert hi["a"] == 255 and lo["a"] == 255
        assert _lum(lo) > _lum(hi) + 0.15, f"low-MAPQ bases should be faded: hi={hi} lo={lo}"

    def test_dv_composite_keeps_strand_stripe(self, viewer_page: Page):
        payload = _payload([_read("fwd", 100), _read("rev", 105, reverse=True)])
        hydrate(viewer_page, payload, "dv-composite")
        # dv-composite rows: 8px tall + 1px gap; stripe is the bottom 2px.
        fwd = sample_pixel_at_genomic(viewer_page, "reads-canvas", 120, 7)
        rev = sample_pixel_at_genomic(viewer_page, "reads-canvas", 120, 9 + 7)
        assert fwd["b"] > fwd["r"] + 60, f"forward composite stripe should be blue: {fwd}"
        assert rev["r"] > rev["b"] + 60, f"reverse composite stripe should be red: {rev}"

    def test_deletion_drawn_as_gap_with_bar(self, viewer_page: Page):
        read = {
            **_read("del", 100, length=40),
            "cigar": "20M5D20M",
            "end_position": 145,
        }
        hydrate(viewer_page, _payload([read]))
        # Mid-row inside the deletion: the dark deletion bar (#1f2937).
        mid = sample_pixel_at_genomic(viewer_page, "reads-canvas", 122, 6)
        assert mid["a"] > 0 and max(mid["r"], mid["g"], mid["b"]) < 80, mid
        # Top edge inside the deletion: cleared (no base painted there).
        top = sample_pixel_at_genomic(viewer_page, "reads-canvas", 122, 1)
        assert top["a"] == 0, top


class TestMismatchLegibility:
    """Pilot-eval regressions: every base colour must stand out on both strands, and
    low-base-quality mismatches must look different from confident ones."""

    # 200 bp in 1024 px → ~5 px/base: wide zoom (ticks), not base tiles.
    WIDE = {"start": 50, "end": 250}
    BODY_Y = 5

    def _wide_read(self, name, reverse, alt, qual=None):
        # Reference at 120 is "A" (ACGTACGTAC repeat from 50); the read carries `alt`.
        r = _read(name, 60, length=120, reverse=reverse)
        r["mismatches"] = [{"pos": 120, "ref": "A", "alt": alt}]
        if qual is not None:
            q = [35] * 120
            q[120 - 60] = qual
            r["qualities"] = q
        return r

    def _tick_and_body(self, page: Page, row: int) -> tuple[dict, dict]:
        y = row * 14 + self.BODY_Y
        return pixel_in_cell(page, 120, y, dx=2), pixel_in_cell(page, 100, y, dx=2)

    @pytest.mark.parametrize("alt", ["A", "C", "G", "T"])
    @pytest.mark.parametrize("reverse", [False, True])
    def test_every_base_visible_on_both_strands(self, viewer_page: Page, alt, reverse):
        hydrate(viewer_page, _payload([self._wide_read("r", reverse, alt)], **self.WIDE))
        tick, body = self._tick_and_body(viewer_page, 0)
        assert _contrast(tick, body) >= 3, (
            f"{alt} tick on {'rev' if reverse else 'fwd'}: {tick} vs {body}"
        )

    def test_low_base_quality_tick_is_hollow(self, viewer_page: Page):
        reads = [
            self._wide_read("hiq", False, "C", qual=35),
            {**self._wide_read("loq", False, "C", qual=8), "position": 60},
        ]
        hydrate(viewer_page, _payload(reads, **self.WIDE))
        rows = viewer_page.evaluate("() => viewer.state.packedRows.map(r => r[0].name)")
        hi_row, lo_row = rows.index("hiq"), rows.index("loq")
        hi_centre = pixel_in_cell(viewer_page, 120, hi_row * 14 + 6, dx=2)
        lo_centre = pixel_in_cell(viewer_page, 120, lo_row * 14 + 6, dx=2)
        lo_cap = pixel_in_cell(viewer_page, 120, lo_row * 14 + 1, dx=2)
        # Confident: solid blue. Low BQ: pale centre, but its border keeps the colour.
        assert lo_centre["b"] > 200 and _lum(lo_centre) > _lum(hi_centre) + 0.3, (
            hi_centre,
            lo_centre,
        )
        assert lo_cap["b"] > lo_cap["r"] + 80, f"low-BQ border should stay blue: {lo_cap}"

    def test_base_zoom_matches_grey_mismatches_coloured(self, viewer_page: Page):
        r = _read("r", 100)
        r["mismatches"] = [{"pos": 120, "ref": "A", "alt": "C"}]
        r["sequence"] = r["sequence"][:20] + "C" + r["sequence"][21:]
        hydrate(viewer_page, _payload([r]))
        match = pixel_in_cell(viewer_page, 110, 1)
        mism = pixel_in_cell(viewer_page, 120, 1)
        assert (
            max(match["r"], match["g"], match["b"]) - min(match["r"], match["g"], match["b"]) < 30
        ), match
        assert mism["b"] > mism["r"] + 100, f"C mismatch tile should be blue: {mism}"
        assert _contrast(mism, match) >= 3
