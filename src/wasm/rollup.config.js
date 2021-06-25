import typescript from "@rollup/plugin-typescript";
import rust from "@wasm-tool/rollup-plugin-rust";

export default {
  input: "index.ts",
  output: {
    file: "dist/sourmash.js",
    format: "umd",
    sourcemap: true,
    name: "sourmash",
  },
  plugins: [rust(), typescript()],
};
