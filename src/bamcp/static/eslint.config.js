// Flat ESLint config for the BAMCP viewer (TypeScript). Minimal but real: JS + TS recommended
// rules. Browser + Node globals come from the `globals` package rather than a hand-rolled list,
// so DOM value-types (HTMLInputElement, localStorage, …) are recognized. typescript-eslint's
// recommended preset already disables core `no-undef` for TS (TypeScript itself checks that);
// the globals still matter for the .js/.ts config files and any non-TS lint.
import js from "@eslint/js";
import globals from "globals";
import tseslint from "typescript-eslint";

export default tseslint.config(
  { ignores: ["dist/**", "node_modules/**"] },
  js.configs.recommended,
  ...tseslint.configs.recommended,
  {
    languageOptions: {
      globals: {
        ...globals.browser,
        ...globals.node,
      },
    },
    rules: {
      "@typescript-eslint/no-explicit-any": "off",
      "@typescript-eslint/no-unused-vars": ["warn", { argsIgnorePattern: "^_" }],
    },
  },
);
