import React from "react";
import HeatMap from "react-heatmap-grid";

export const Heatmap = () => {
  const xLabels = [
    "HP:0020225", // Breast abscess
    "HP:6000671", // Breast myxoma
    "HP:6000894", // Breast adenoma
    "HP:6001046", // Breast cyst
    "HP:0000769", // Abnormality of the breast
    "HP:4000122", // History of exclusive breast feeding
    "HP:6000895", // Breast apocrine adenoma
    "HP:0010313", // Breast hypertrophy
    "HP:6000102", // Breast intraductal papilloma
    "HP:0030901", // Pruritis on breast
    "HP:0000769", // Breast aplasia
    "HP:6000102", // Breast mass
    "HP:0010313", // Mastalgia
    "HP:6000895", // Architectural distortion of breast
    "HP:0000769", // Corneal transplant history
    "HP:0010313", // Unilateral breast hypoplasia
    "HP:6001046", // Hypoechoic breast mass
  ];

  const yLabels = [
    "OR4F5", // Uncertain significance
    "CFTR",  // Pathogenic
    "BRCA1", // Likely pathogenic
    "DMD",   // Pathogenic
    "FGFR3", // Benign
  ];

  const dataMatrix = [
    // OR4F5 - Uncertain significance
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2],

    // CFTR - Pathogenic
    [3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 3, 3],

    // BRCA1 - Likely pathogenic
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 3, 3],

    // DMD - Pathogenic
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],

    // FGFR3 - Benign
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1],
  ];

  return (
    <div style={{ fontSize: "12px", width: "100%", overflowX: "auto" }}>
      <h2>Patient Phenotypic and Genomic Heatmap</h2>
      <div style={{ padding: "20px 0 0 100px" }}>
        <HeatMap
          xLabels={xLabels}
          yLabels={yLabels}
          data={dataMatrix}
          xLabelsLocation="top"
          cellStyle={(background, value, min, max, data, x, y) => ({
            background: ["#ffffff", "#ffeda0", "#feb24c", "#f03b20"][value] || "#ffffff",
            fontSize: "11px",
            border: "1px solid #ccc",
            padding: "4px",
            height: "30px",
          })}
          xLabelsStyle={() => ({
            color: "black",
            fontSize: "11px",
            transform: "rotate(-45deg)",
            transformOrigin: "left",
            textAlign: "left",
            paddingBottom: "10px",
            paddingLeft: "5px",
            height: "100px",
            width: "100px",
          })}
          yLabelsStyle={() => ({
            fontSize: "11px",
            textAlign: "right",
            paddingRight: "10px",
          })}
          cellRender={(value) => (value !== 0 ? value : "")}
        />
      </div>
    </div>
  );
};
