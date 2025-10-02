# 🦠 Interactive SEIRS Model

This interactive simulator demonstrates how infectious diseases spread through populations using a deterministic, mean-field SEIRS compartmental model. 
It features real-time parameter adjustments and elegant visualizations to help understand epidemic dynamics, including vaccination, immunity waning, 
and disease-induced mortality.

<br>

<div align="center">

**[🌐 Live Demo](https://denniemok.github.io/seirs-demo)** • **[📥 Download](https://github.com/denniemok/seirs-demo/archive/refs/heads/main.zip)** • **[📖 Documentation](#-technical-details)**

</div>

<br>

## ✨ Features

<table>
<tr>
<td width=500px>

### 🎮 Interactive Controls
- Real-time parameter adjustments
- Intuitive slider controls
- Instant visualization updates
- Reset to defaults anytime
<br>

</td>
<td width=500px>

### 📊 Flexible Visualization
- Linear & logarithmic Y-axis
- Adjustable time horizons (3000 days)
- Color-coded compartments (S, E, I, R)
- Responsive plot design
<br>

</td>
</tr>
<tr>
<td >

### 🔬 Epidemiological Parameters
- Basic reproduction number (R₀)
- Incubation & infectious periods
- Immunity duration
- Vaccination rates
- Disease-induced mortality
<br>

</td>
<td >

### 🎓 Educational Content
- Built-in mathematical equations
- Parameter explanations
- Helpful tooltips
- Clean, modern interface
- Zero dependencies required
<br>

</td>
</tr>
</table>

<br>

## 🚀 Getting Started

### Prerequisites

- A modern web browser (Chrome, Firefox, Safari, or Edge)
- No installation or build process required!

### Quick Start

1. #### Get the code

    ```bash
    git clone https://github.com/denniemok/seirs-demo.git
    ```
    
    Or [download ZIP](https://github.com/denniemok/seirs-demo/archive/refs/heads/main.zip) directly.

2. #### Open in browser

    Double-click `index.html` to open directly in your browser, or use a local web server for full functionality:
    
    ```bash
    cd seirs-demo
    
    # Python 3
    python -m http.server 8000
    
    # Node.js
    npx http-server
    ```

    Then open **`http://localhost:8000`** in your browser.

</details>

<br>

## 📐 Technical Details

### Framework

The model divides the population into **four compartments**:

| Compartment | Description |
|-------------|-------------|
| **S** (Susceptible) | 🟦 Individuals who can contract the disease |
| **E** (Exposed) | 🟩 Individuals who are infected but not yet infectious |
| **I** (Infectious) | 🟧 Individuals who can transmit the disease |
| **R** (Recovered) | 🟨 Individuals who have immunity (temporary or permanent) |

<br>

### Mathematics

The dynamics are governed by these differential equations:

```
dS/dt = -βSI + ωR - μS + μ(1-p)
dE/dt = βSI - σE - μE
dI/dt = σE - γI - (μ+α)I
dR/dt = γI - ωR - μR + μp
```

| Symbol | Description |
|--------|-------------|
| **β** | Transmission rate |
| **σ** | Progression rate (1 / incubation period) |
| **γ** | Recovery rate (1 / infectious period) |
| **ω** | Immunity waning rate (1 / immunity duration) |
| **μ** | Natural mortality rate (1 / life expectancy) |
| **α** | Disease-induced mortality rate (1 / infection-to-death period) |
| **p** | Vaccination rate |

<br>

### Implementation

| Feature | Technology |
|---------|------------|
| Numerical solver | 4th-order Runge-Kutta (RK4) |
| Data structures | `Float64Array` for performance |
| Visualization | D3.js v3.5.17 (SVG) |
| Documentation | JSDoc + inline comments |

<br>

## 📁 Project Structure

```
seirs-demo/
├── 📄 index.html          # Main HTML page with structure and equations
├── 📜 seirs.js            # Core SEIRS model implementation and plotting logic
├── ⚙️  params.js           # Parameter definitions and UI control values
├── 🎨 seirs.css           # Styling and responsive design
└── 📊 d3.min.js           # D3.js library for visualization
```

<br>

## 🎨 Customisation

### Modify Parameters

Edit `params.js` to change parameter ranges and defaults:

```javascript
// Example: Extend R₀ range from 1-5 to 1-10
R0: generateParams(1, 10, 0.1, 3.0)
```

### Change Styling

Modify `seirs.css` to customize the appearance. CSS variables make theming easy:

```css
:root {
    --color-susceptible: #2c3e50;  /* Dark blue */
    --color-exposed: #27ae60;      /* Green */
    --color-infectious: #e67e22;   /* Orange */
    --color-recovered: #3498db;    /* Light blue */
}
```

<br>

## 🙏 Acknowledgments

This project is built upon and inspired by:

- 📚 [seir-demo](https://github.com/robmoss/seir-demo) by Rob Moss
- 📚 [posepi2](https://github.com/martinkrz/posepi2) by Martin Krzywinski
- 📊 [D3.js](https://d3js.org/) by Mike Bostock
- 📖 Bjørnstad, O., Shea, K., Krzywinski, M. & Altman, N. [_The SEIRS model for infectious disease dynamics_](http://www.nature.com/articles/s41592-020-0856-2). Nature Methods **17**:557–558 (2020)

<br>
<hr>

This model is intended for educational purposes to understand epidemiological dynamics. For public health decisions, please consult epidemiological experts and use validated, peer-reviewed models.

This project, including D3.js v3.5.17, is distributed under the [**BSD 3-Clause License**](LICENSE).

<br><br><br>

<div align="center">

Made with ❤️ for epidemiology education

**[⬆ Back to Top](#-interactive-seirs-model)**

</div>
