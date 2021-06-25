export * from "./target/wasm-pack/sourmash-wasm/index";

type Exports = typeof import("./target/wasm-pack/sourmash-wasm/index");
declare const init: () => Promise<Exports>;
export default init;
