import React, { useState, useEffect } from "react";
import { Flame, ShieldAlert, AlertCircle } from "lucide-react";

export interface ReactiveInteractionDetail {
  row: number;
  col: number;
  rowLabel: string;
  colLabel: string;
  score: number;
  percentage: string;
  mechanism: string;
  severity: "Critical" | "High" | "Moderate" | "Low";
  description: string;
}

export interface HeatmapResponse {
  success: boolean;
  compound1Name: string;
  compound2Name: string;
  rowLabels: string[];
  colLabels: string[];
  matrix: number[][];
  details: ReactiveInteractionDetail[];
  topInteractions: ReactiveInteractionDetail[];
  heatmapBase64: string;
  summary: string;
  error?: string;
}

interface InteractionHeatmapProps {
  compounds: Array<{
    name: string;
    smiles?: string;
    features?: string[];
    interactionSites?: string[];
  }>;
}

export const InteractionHeatmap: React.FC<InteractionHeatmapProps> = ({ compounds }) => {
  const [loading, setLoading] = useState<boolean>(false);
  const [data, setData] = useState<HeatmapResponse | null>(null);
  const [error, setError] = useState<string | null>(null);

  const primaryCompound = compounds[0];
  const secondaryCompounds = compounds.slice(1);
  const hasSecondaries = secondaryCompounds.length > 0;

  // Target pair: Primary vs Secondary (or Intramolecular if solo)
  const currentPair = hasSecondaries
    ? [primaryCompound, secondaryCompounds[0]]
    : [primaryCompound];

  useEffect(() => {
    if (!compounds || compounds.length === 0) return;
    let isSubscribed = true;
    setLoading(true);
    setError(null);

    fetch("/api/interaction-heatmap", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        compounds: currentPair,
        cmap: "warmcool",
        title: hasSecondaries
          ? `Stress Degradation & Incompatibility: ${currentPair[0].name} & ${currentPair[1].name}`
          : `Stress Degradation Profile: ${currentPair[0].name}`,
      }),
    })
      .then(async (res) => {
        if (!res.ok) {
          const errJson = await res.json().catch(() => ({}));
          throw new Error(errJson.error || `Failed with status ${res.status}`);
        }
        return res.json();
      })
      .then((json: HeatmapResponse) => {
        if (!isSubscribed) return;
        if (!json.success) {
          throw new Error(json.error || "Heatmap calculation encountered an issue.");
        }
        setData(json);
      })
      .catch((err: any) => {
        if (!isSubscribed) return;
        console.error("Heatmap fetch error:", err);
        setError(err?.message || "Failed to generate interaction heatmap.");
      })
      .finally(() => {
        if (isSubscribed) setLoading(false);
      });

    return () => {
      isSubscribed = false;
    };
  }, [compounds]);

  return (
    <section className="bg-white border border-[#E2E8F0] rounded-2xl p-6 shadow-xs space-y-6">
      {/* Header - Completely static */}
      <div>
        <div className="flex items-center gap-2">
          <div className="w-8 h-8 rounded-lg bg-[#EFF6FF] flex items-center justify-center text-[#2563EB]">
            <Flame className="w-4 h-4" />
          </div>
          <h3 className="font-serif text-xl font-bold text-[#0F172A]">
            Stress Degradation &amp; Incompatibility Heatmap
          </h3>
        </div>
        <p className="text-xs sm:text-sm text-[#64748B] mt-1">
          Quantitative stress matrix modeling reactive center vulnerability across Acidic, Basic, Hydrolysis, Photolysis, Thermal, and Oxidative conditions using the WarmCool spectrum.
        </p>
      </div>

      {/* Error View */}
      {error && (
        <div className="p-4 bg-amber-50 border border-amber-200 rounded-xl flex items-start gap-3 text-amber-900 text-xs">
          <AlertCircle className="w-4 h-4 text-amber-600 mt-0.5 shrink-0" />
          <div>
            <div className="font-bold mb-0.5">Heatmap Rendering Notice</div>
            <div>{error}</div>
          </div>
        </div>
      )}

      {/* Loading State */}
      {loading && (
        <div className="py-16 flex flex-col items-center justify-center text-center space-y-3 bg-[#F8FAFC] rounded-xl border border-dashed border-[#CBD5E1]">
          <div className="w-8 h-8 border-3 border-[#EFF6FF] border-t-[#2563EB] rounded-full animate-spin" />
          <div className="text-sm font-semibold text-[#0F172A]">
            Generating WarmCool Stress Degradation Matrix...
          </div>
          <div className="text-xs text-[#64748B] max-w-sm">
            Evaluating reactive center vulnerability across Acidic, Basic, Hydrolysis, Photolysis, Thermal, and Oxidative stress conditions.
          </div>
        </div>
      )}

      {/* Main Heatmap Visual Display - Static (Reactive Centers Cross-Affinity Matrix section removed) */}
      {!loading && data?.heatmapBase64 && (
        <div className="space-y-6">
          <div className="bg-[#F8FAFC] p-4 rounded-xl border border-[#E2E8F0] flex flex-col items-center justify-center">
            <img
              src={data.heatmapBase64}
              alt="WarmCool Stress Degradation Heatmap"
              className="max-h-[500px] w-auto object-contain rounded-lg shadow-xs"
            />
            <div className="mt-3 text-xs text-[#64748B] text-center italic max-w-2xl">
              {data.summary}
            </div>
          </div>

          {/* Static Top Reactive Couplings & Degradation Pathways Ranking */}
          {data.topInteractions && data.topInteractions.length > 0 && (
            <div className="space-y-3">
              <div className="flex items-center justify-between">
                <div className="text-xs font-bold text-[#64748B] uppercase tracking-wider flex items-center gap-1.5">
                  <ShieldAlert className="w-4 h-4 text-[#DC2626]" />
                  Critical Degradation &amp; Incompatibility Conditions (Ranked by Reaction Potential)
                </div>
                <span className="text-xs text-[#94A3B8]">Top {data.topInteractions.length} Exposure Points</span>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                {data.topInteractions.map((item, idx) => {
                  const isCrit = item.severity === "Critical";
                  const isHigh = item.severity === "High";

                  const cardBorder = isCrit
                    ? "bg-[#FEF2F2] border-[#FECACA] text-[#7F1D1D]"
                    : isHigh
                    ? "bg-[#FFF7ED] border-[#FED7AA] text-[#7C2D12]"
                    : "bg-[#EFF6FF] border-[#BFDBFE] text-[#1E3A8A]";

                  const badgeBg = isCrit
                    ? "bg-[#DC2626] text-white"
                    : isHigh
                    ? "bg-[#EA580C] text-white"
                    : "bg-[#2563EB] text-white";

                  return (
                    <div
                      key={`top-int-${idx}`}
                      className={`p-4 rounded-xl border ${cardBorder}`}
                    >
                      <div className="flex items-center justify-between mb-2">
                        <div className="font-bold text-xs">
                          {item.colLabel} &bull; <span className="underline decoration-[#94A3B8]">{item.rowLabel} Stress</span>
                        </div>
                        <span className={`font-mono font-bold text-xs px-2 py-0.5 rounded-full ${badgeBg}`}>
                          {item.severity}
                        </span>
                      </div>
                      <div className="text-xs font-semibold text-[#1D4ED8] mb-1">
                        {item.mechanism}
                      </div>
                      <p className="text-xs text-[#475569] leading-relaxed">
                        {item.description}
                      </p>
                    </div>
                  );
                })}
              </div>
            </div>
          )}
        </div>
      )}
    </section>
  );
};
