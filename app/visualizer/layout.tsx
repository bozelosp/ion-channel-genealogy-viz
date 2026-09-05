import { preload } from "react-dom";

export default function VisualizerLayout({ children }: Readonly<{ children: React.ReactNode }>) {
  preload("/data/network_data_initial.json", { as: "fetch", crossOrigin: "anonymous" });
  return children;
}
