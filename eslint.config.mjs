import { defineConfig, globalIgnores } from "eslint/config";
import nextVitals from "eslint-config-next/core-web-vitals";
import nextTypescript from "eslint-config-next/typescript";

const eslintConfig = defineConfig([
  ...nextVitals,
  ...nextTypescript,
  {
    rules: {
      // Preserve the established effect-driven synchronization semantics. These
      // call sites require behavioral refactors, not a security-upgrade rewrite.
      "react-hooks/set-state-in-effect": "off",
    },
  },
  globalIgnores([
    "node_modules/**",
    ".next/**",
    "out/**",
    "build/**",
    "public/**",
    "app/omnimodels/static/generated/**",
    "**/node_modules/**",
    "next-env.d.ts",
  ]),
]);

export default eslintConfig;
