import { spawn } from "child_process";
import path from "path";

export interface SeabornHeatmapInput {
  matrix: number[][];
  rowLabels: string[];
  colLabels: string[];
  title?: string;
  cmap?: string;
}

/**
 * Invokes python3 scripts/plot_seaborn_heatmap.py to produce a high-resolution,
 * publication-grade Seaborn heatmap PNG (base64 data URL).
 */
export async function renderSeabornHeatmap(data: SeabornHeatmapInput): Promise<string | null> {
  return new Promise((resolve) => {
    try {
      const scriptPath = path.join(process.cwd(), "scripts", "plot_seaborn_heatmap.py");
      const child = spawn("python3", [scriptPath]);
      let stdout = "";
      let stderr = "";

      child.stdout.on("data", (chunk) => {
        stdout += chunk;
      });

      child.stderr.on("data", (chunk) => {
        stderr += chunk;
      });

      child.on("close", (code) => {
        if (code === 0 && stdout) {
          try {
            const parsed = JSON.parse(stdout.trim());
            if (parsed.success && parsed.image) {
              return resolve(parsed.image);
            }
          } catch (e) {
            console.warn("Could not parse Python Seaborn output:", e);
          }
        }
        if (stderr) {
          console.warn("Python Seaborn warning/stderr:", stderr);
        }
        resolve(null);
      });

      child.on("error", (err) => {
        console.warn("Python process spawn error:", err.message);
        resolve(null);
      });

      child.stdin.write(JSON.stringify(data));
      child.stdin.end();
    } catch (err) {
      console.warn("Failed to launch Seaborn script:", err);
      resolve(null);
    }
  });
}
