import { createRoot } from "react-dom/client";
import App from "./App";
import "./styles.css";

// no StrictMode: the 3D viewer owns a WebGL canvas that must be created exactly once
createRoot(document.getElementById("root")!).render(<App />);
