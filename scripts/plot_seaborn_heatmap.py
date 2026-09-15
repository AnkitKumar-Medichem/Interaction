import sys
import json
import base64
import io

def main():
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
        import numpy as np
        import pandas as pd
    except ImportError as e:
        print(json.dumps({"success": False, "error": f"Missing python library: {str(e)}"}))
        sys.exit(1)

    try:
        if len(sys.argv) > 1:
            with open(sys.argv[1], 'r') as f:
                data = json.load(f)
        else:
            data = json.load(sys.stdin)

        matrix = np.array(data.get('matrix', []))
        row_labels = data.get('rowLabels', [])
        col_labels = data.get('colLabels', [])
        title = data.get('title', 'Stress Degradation & Incompatibility Heatmap')
        cmap_name = data.get('cmap', 'coolwarm')

        if matrix.size == 0 or len(row_labels) == 0 or len(col_labels) == 0:
            print(json.dumps({"success": False, "error": "Empty matrix data provided"}))
            sys.exit(1)

        # Truncate long row labels for clean display
        display_row_labels = [lbl if len(lbl) <= 38 else lbl[:35] + "..." for lbl in row_labels]

        df = pd.DataFrame(matrix, index=display_row_labels, columns=col_labels)

        # Dimensions
        n_rows = len(display_row_labels)
        n_cols = len(col_labels)
        fig_width = max(9.0, n_cols * 1.5)
        fig_height = max(5.2, n_rows * 0.85 + 1.8)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=150)
        fig.patch.set_facecolor('#FFFFFF')
        ax.set_facecolor('#F8FAFC')

        # Annotation labels formatted as percentages
        annot_matrix = np.vectorize(lambda x: f"{int(round(float(x) * 100))}%")(matrix)

        # Diverging colormap
        cmap = 'vlag' if cmap_name in ['vlag', 'warmcool'] else 'coolwarm'

        sns.heatmap(
            df,
            annot=annot_matrix,
            fmt="",
            cmap=cmap,
            vmin=0.0,
            vmax=1.0,
            cbar_kws={'label': 'Degradation / Incompatibility Potential', 'shrink': 0.85},
            linewidths=2.0,
            linecolor='#FFFFFF',
            square=False,
            ax=ax,
            annot_kws={'fontsize': 10, 'fontweight': 'bold'}
        )

        ax.set_title(title, fontsize=13, fontweight='bold', pad=18, color='#0F172A')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=20, ha='right', fontsize=9.5, fontweight='600', color='#334155')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9.5, fontweight='600', color='#334155')

        plt.tight_layout()

        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight', dpi=180)
        plt.close(fig)
        buf.seek(0)
        b64_str = base64.b64encode(buf.read()).decode('utf-8')
        
        print(json.dumps({
            "success": True, 
            "image": f"data:image/png;base64,{b64_str}"
        }))
    except Exception as e:
        print(json.dumps({"success": False, "error": str(e)}))
        sys.exit(1)

if __name__ == '__main__':
    main()
