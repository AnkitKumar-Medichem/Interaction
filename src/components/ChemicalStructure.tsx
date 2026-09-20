import React, { useEffect, useState } from 'react';
import { getMoleculeSvg, sanitizeSmiles } from '../lib/rdkit';

interface ChemicalStructureProps {
  smiles: string;
  width?: number;
  height?: number;
  className?: string;
  altText?: string;
}

export const ChemicalStructure: React.FC<ChemicalStructureProps> = ({ 
  smiles, 
  width = 200, 
  height = 200,
  className = "",
  altText = "Chemical Structure"
}) => {
  const [svgContent, setSvgContent] = useState<string | null>(null);
  const [fallbackImgUrl, setFallbackImgUrl] = useState<string | null>(null);
  const [attemptStage, setAttemptStage] = useState<number>(0);
  const [isLoading, setIsLoading] = useState<boolean>(true);
  const [error, setError] = useState<boolean>(false);

  useEffect(() => {
    let mounted = true;
    const clean = sanitizeSmiles(smiles);
    
    if (!clean) {
      setError(true);
      setIsLoading(false);
      setSvgContent(null);
      setFallbackImgUrl(null);
      return;
    }

    setError(false);
    setIsLoading(true);
    setSvgContent(null);
    setFallbackImgUrl(null);
    setAttemptStage(0);

    // Stage 1: Try in-browser RDKit WASM (fastest, client-side vector SVG)
    getMoleculeSvg(clean, width, height)
      .then((svg) => {
        if (!mounted) return;
        if (svg) {
          setSvgContent(svg);
          setFallbackImgUrl(null);
          setIsLoading(false);
          setError(false);
        } else {
          // Stage 2: In-browser RDKit couldn't generate SVG, try server-side depiction API
          tryServerDepiction(clean, mounted);
        }
      })
      .catch(() => {
        if (!mounted) return;
        tryServerDepiction(clean, mounted);
      });

    return () => { mounted = false; };
  }, [smiles, width, height]);

  const tryServerDepiction = (clean: string, mounted: boolean) => {
    fetch(`/api/structure?smiles=${encodeURIComponent(clean)}&w=${width}&h=${height}`)
      .then(async (res) => {
        if (!mounted) return;
        if (res.ok) {
          const contentType = res.headers.get("content-type") || "";
          if (contentType.includes("svg")) {
            const svgText = await res.text();
            if (svgText && svgText.includes("<svg")) {
              setSvgContent(svgText);
              setFallbackImgUrl(null);
              setIsLoading(false);
              setError(false);
              return;
            }
          } else {
            // Received image bytes (e.g. from Cactus proxy)
            const blob = await res.blob();
            const objectUrl = URL.createObjectURL(blob);
            setFallbackImgUrl(objectUrl);
            setIsLoading(false);
            setError(false);
            return;
          }
        }
        // Stage 3: Try direct Cactus NCI depiction
        fallbackToExternalCactus(clean, mounted);
      })
      .catch(() => {
        if (!mounted) return;
        fallbackToExternalCactus(clean, mounted);
      });
  };

  const fallbackToExternalCactus = (clean: string, mounted: boolean) => {
    if (!mounted) return;
    setAttemptStage(1);
    setSvgContent(null);
    setFallbackImgUrl(`https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(clean)}/image`);
    setIsLoading(false);
  };

  const handleImgError = () => {
    const clean = sanitizeSmiles(smiles);
    if (attemptStage === 1) {
      // Cactus failed, try PubChem REST depiction
      setAttemptStage(2);
      setFallbackImgUrl(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(clean)}/PNG?record_type=2d&image_size=${width}x${height}`);
    } else {
      // All image depiction attempts exhausted
      setError(true);
      setFallbackImgUrl(null);
      setIsLoading(false);
    }
  };

  return (
    <div 
      className={`relative flex items-center justify-center bg-white rounded-lg overflow-hidden border border-slate-100 ${className}`}
      style={{ width, height }}
    >
      {isLoading && (
        <div className="absolute inset-0 flex flex-col items-center justify-center bg-slate-50/70 animate-pulse">
          <div className="w-8 h-8 rounded-full border-2 border-slate-200 border-t-blue-500 animate-spin mb-1.5" />
          <span className="text-[10px] text-slate-400 font-mono">Rendering...</span>
        </div>
      )}

      {error && !isLoading && (
        <div className="absolute inset-0 flex flex-col items-center justify-center bg-slate-50 text-slate-400 p-2 text-center select-none">
          <svg className="w-6 h-6 mb-1 text-slate-300" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
          </svg>
          <span className="text-[10px] font-medium text-slate-500">2D Diagram</span>
          <span className="text-[9px] text-slate-400 truncate max-w-[90%] font-mono mt-0.5" title={smiles}>
            {sanitizeSmiles(smiles).slice(0, 18)}...
          </span>
        </div>
      )}

      {!error && svgContent && (
        <div 
          className="w-full h-full flex items-center justify-center p-1.5 [&>svg]:max-w-full [&>svg]:max-h-full [&>svg]:w-auto [&>svg]:h-auto transition-opacity duration-200"
          dangerouslySetInnerHTML={{ __html: svgContent }} 
        />
      )}

      {!error && !svgContent && fallbackImgUrl && (
        <img
          src={fallbackImgUrl}
          alt={altText}
          className="max-w-full max-h-full object-contain p-2 transition-opacity duration-200"
          onError={handleImgError}
          onLoad={() => setIsLoading(false)}
          referrerPolicy="no-referrer"
        />
      )}
    </div>
  );
};
