// components/HeatmapComponent.js

import React from "react";
import HeatMap from "react-heatmap-grid";

export const Heatmap = () => {
  const xLabels = [
    "Breast abscess",
    "Breast myxoma",
    "Breast adenoma",
    "Breast cyst",
    "Abnormality of the breast",
    "History of exclusive breast feeding",
    "Breast apocrine adenoma",
    "Breast hypertrophy",
    "Breast intraductal papilloma",
    "Pruritis on breast",
    "Breast aplasia",
    "Breast mass",
    "Mastalgia",
    "Architectural distortion of breast",
    "Corneal transplant history",
    "Unilateral breast hypoplasia",
    "Hypoechoic breast mass",
  ];

  const yLabels = [
    "OR4F5 (missense_variant)",
    "CFTR (frameshift_variant)",
    "BRCA1 (stop_gained)",
    "DMD (splice_donor_variant)",
    "FGFR3 (missense_variant)",
  ];

  const dataMatrix = [
    // OR4F5
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],

    // CFTR
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],

    // BRCA1
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 3, 3],

    // DMD
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],

    // FGFR3
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
  ];

  return (
    <div style={{ fontSize: "12px", overflowX: "auto" }}>
      <h2>Patient Phenotypic and Genomic Heatmap</h2>
      <HeatMap
        xLabels={xLabels}
        yLabels={yLabels}
        data={dataMatrix}
        squares
        height={45}
        cellStyle={(background, value, min, max, data, x, y) => {
          const colors = ["#ffffff", "#ffeda0", "#feb24c", "#f03b20"];
          return {
            backgroundColor: colors[value] || "#ffffff",
            fontSize: "12px",
            border: "1px solid #ccc",
          };
        }}
        xLabelsStyle={(index) => ({
          color: "black",
          fontSize: "12px",
          transform: "rotate(-90deg)",
          textAlign: "left",
          width: "150px",
          display: "block",
          overflow: "hidden",
          whiteSpace: "nowrap",
        })}
        yLabelsStyle={() => ({
          fontSize: "12px",
          textAlign: "right",
          width: "250px",
        })}
        cellRender={(value) => (value !== 0 ? value : "")}
      />
    </div>
  );
};
