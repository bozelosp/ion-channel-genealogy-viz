import { defineConfig, globalIgnores } from "eslint/config";
import js from "@eslint/js";
import jsxA11y from "eslint-plugin-jsx-a11y";
import react from "eslint-plugin-react";
import reactHooks from "eslint-plugin-react-hooks";
import tseslint from "typescript-eslint";

const eslintConfig = defineConfig([
  js.configs.recommended,
  ...tseslint.configs.recommended,
  react.configs.flat.recommended,
  react.configs.flat["jsx-runtime"],
  jsxA11y.flatConfigs.recommended,
  reactHooks.configs.flat.recommended,
  {
    settings: { react: { version: "detect" } },
    rules: {
      // Existing browser fallbacks intentionally ignore a few optional failures.
      "no-empty": "off",
      // Context-menu propagation, dialog key trapping, and a scrollable <pre>
      // are deliberate non-control interactions.
      "jsx-a11y/click-events-have-key-events": "off",
      "jsx-a11y/no-static-element-interactions": "off",
      "jsx-a11y/no-noninteractive-element-interactions": "off",
      "jsx-a11y/no-noninteractive-tabindex": "off",
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
