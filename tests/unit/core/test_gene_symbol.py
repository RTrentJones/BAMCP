"""Gene-symbol validation at the search_gene boundary and in the NCBI client."""

import pytest

from bamcp.clients.genes import GeneClient
from bamcp.config import BAMCPConfig
from bamcp.core.tools import handle_search_gene
from bamcp.core.validation import validate_gene_symbol

_BAD = [
    "BRCA1[Gene Name] OR 1[uid]",
    "TP53 OR BRCA1",
    'TP53"',
    "TP53)",
    "-TP53",
    ".x",
    "a" * 41,
    "TP53\n",
    "TP53 ",
    "BRCA1#",
    "",
]


@pytest.mark.unit
@pytest.mark.parametrize(
    "symbol", ["BRCA1", "TP53", "HLA-A", "C4orf3", "MT-ND1", "A1BG.AS1", "a" * 40]
)
def test_valid_symbols(symbol):
    validate_gene_symbol(symbol)


@pytest.mark.unit
@pytest.mark.parametrize("symbol", _BAD)
def test_invalid_symbols(symbol):
    with pytest.raises(ValueError, match="Invalid gene symbol"):
        validate_gene_symbol(symbol)


@pytest.mark.unit
async def test_handler_rejects_before_any_lookup(httpx_mock):
    with pytest.raises(ValueError, match="Invalid gene symbol"):
        await handle_search_gene({"symbol": "BRCA1[Gene Name] OR x"}, BAMCPConfig())
    assert httpx_mock.get_requests() == []


@pytest.mark.unit
async def test_client_rejects_defensively(httpx_mock):
    client = GeneClient()
    with pytest.raises(ValueError, match="Invalid gene symbol"):
        await client.search("TP53 OR BRCA1")
    assert httpx_mock.get_requests() == []
