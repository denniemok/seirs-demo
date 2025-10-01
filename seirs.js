// SEIRS model implementation

// ============================================================================
// Constants
// ============================================================================

const CONSTANTS = {
    DAYS_PER_YEAR: 365,               // Conversion factor for yearly rates to daily rates
    PERCENTAGE_SCALE: 100,            // Scale factor for displaying proportions as percentages
    RK4: {
        WEIGHT_MIDDLE: 2,             // Weight for k2 and k3 in RK4 formula
        WEIGHT_DIVISOR: 6,            // Divisor in RK4 weighted average
        HALF_STEP: 0.5                // Half step size for intermediate RK4 calculations
    }
};

// ============================================================================
// Validation Utilities
// ============================================================================

/**
 * Validates input parameters for the SEIR model
 */
class ParameterValidator {
    /**
     * Validates that a value is within a specified range
     * @param {number} value - The value to validate
     * @param {number} min - Minimum allowed value
     * @param {number} max - Maximum allowed value
     * @param {string} name - Parameter name for error messages
     */
    static validateRange(value, min, max, name) {
        if (value < min || value > max) {
            throw new Error(`${name} must be between ${min} and ${max}`);
        }
    }

    /**
     * Validates that a value is positive (> 0)
     */
    static validatePositive(value, name) {
        if (value <= 0) {
            throw new Error(`${name} must be positive`);
        }
    }

    /**
     * Validates that a value is non-negative (>= 0)
     */
    static validateNonNegative(value, name) {
        if (value < 0) {
            throw new Error(`${name} must be non-negative`);
        }
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Clamps a value between min and max bounds
 * Used to ensure compartment values stay in valid range [0,1]
 */
const clamp = (value, min = 0, max = 1) => Math.max(min, Math.min(max, value));

/**
 * Formats numerical array data for plotting
 * Converts proportions to percentages and creates {x, y} coordinate pairs
 */
const formatDataForPlot = (array, scale = CONSTANTS.PERCENTAGE_SCALE) => {
    return Array.from(array, (value, index) => ({
        x: index,
        y: scale * value
    }));
};

// ============================================================================
// SEIR Model Core
// ============================================================================

/**
 * Calculates the rate parameters for the SEIR model
 * 
 * Model parameters:
 * - beta: Transmission rate (contact rate × probability of transmission)
 * - sigma: Rate of progression from exposed to infectious (1/latent_period)
 * - gamma: Recovery rate (1/infectious_period)
 * - alpha: Disease-induced death rate (1/death_onset)
 * - mu: Natural death/birth rate (1/life_expectancy)
 * - omega: Immunity waning rate (1/immunity_duration)
 */
class SEIRParameters {
    constructor(config) {
        const { R0, infectious_period, latent_period, death_onset, immunity_duration, life_expectancy } = config;
        
        // Disease-induced death rate (0 if death_onset is infinite)
        this.alpha = death_onset > 0 ? 1 / death_onset : 0;
        
        // Recovery rate (1/days to recover)
        this.gamma = 1 / infectious_period;
        
        // Immunity waning rate (convert years to days)
        this.omega = 1 / (CONSTANTS.DAYS_PER_YEAR * immunity_duration);
        
        // Natural death/birth rate (convert years to days)
        this.mu = 1 / (CONSTANTS.DAYS_PER_YEAR * life_expectancy);
        
        // Rate of progression from exposed to infectious
        this.sigma = 1 / latent_period;
        
        // Transmission rate calculated from R0 and other parameters
        // R0 = beta * sigma / ((sigma + mu) * (gamma + mu + alpha))
        this.beta = R0 * (this.gamma + this.mu + this.alpha) * (this.sigma + this.mu) / this.sigma;
    }
}

/**
 * Calculates the rate of change (derivatives) for each SEIR compartment
 * 
 * Differential equations:
 * dS/dt = -beta*S*I + omega*R - mu*S + mu*(1-p)
 * dE/dt = beta*S*I - sigma*E - mu*E
 * dI/dt = sigma*E - gamma*I - (mu+alpha)*I
 * dR/dt = gamma*I - omega*R - mu*R + mu*p
 * 
 * where p is the vaccination rate
 */
class SEIRTransitions {
    constructor(params, vaccination_rate) {
        this.params = params;
        this.vaccination_rate = vaccination_rate;
    }

    /**
     * Calculates derivatives for all compartments at current state
     * @param {number} s - Susceptible proportion
     * @param {number} e - Exposed proportion
     * @param {number} i - Infectious proportion
     * @param {number} r - Recovered proportion
     * @returns {Object} Derivatives {ds, de, di, dr}
     */
    calculateDerivatives(s, e, i, r) {
        const { beta, sigma, gamma, omega, mu, alpha } = this.params;
        const p = this.vaccination_rate;

        // Calculate transition rates between compartments
        const s_to_e = beta * s * i;      // Force of infection
        const e_to_i = sigma * e;          // Progression to infectious
        const i_to_r = gamma * i;          // Recovery
        const r_to_s = omega * r;          // Immunity waning

        return {
            // Susceptible: lose to infection, gain from waning immunity and births
            ds: -s_to_e + r_to_s - mu * s + mu * (1 - p),
            
            // Exposed: gain from infection, lose to progression and death
            de: s_to_e - e_to_i - mu * e,
            
            // Infectious: gain from exposed, lose to recovery and death
            di: e_to_i - i_to_r - (mu + alpha) * i,
            
            // Recovered: gain from recovery and vaccinated births, lose to waning and death
            dr: i_to_r - r_to_s - mu * r + mu * p
        };
    }
}

/**
 * Runge-Kutta 4th order numerical integrator
 * 
 * RK4 is a high-accuracy method for solving ODEs:
 * x(t+h) = x(t) + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
 * 
 * where:
 * k1 = f(x, t)
 * k2 = f(x + 0.5*h*k1, t + 0.5*h)
 * k3 = f(x + 0.5*h*k2, t + 0.5*h)
 * k4 = f(x + h*k3, t + h)
 */
class RK4Integrator {
    constructor(transitions) {
        this.transitions = transitions;
    }

    /**
     * Performs one RK4 integration step
     * @param {Array} state - Current state [s, e, i, r]
     * @param {number} h - Time step size
     * @returns {Array} New state [s_new, e_new, i_new, r_new]
     */
    step(state, h) {
        const [s, e, i, r] = state;
        const { HALF_STEP, WEIGHT_MIDDLE, WEIGHT_DIVISOR } = CONSTANTS.RK4;

        // k1: slope at beginning of interval
        const k1 = this.transitions.calculateDerivatives(s, e, i, r);
        
        // k2: slope at midpoint using k1
        const k2 = this.transitions.calculateDerivatives(
            s + HALF_STEP * h * k1.ds,
            e + HALF_STEP * h * k1.de,
            i + HALF_STEP * h * k1.di,
            r + HALF_STEP * h * k1.dr
        );

        // k3: slope at midpoint using k2
        const k3 = this.transitions.calculateDerivatives(
            s + HALF_STEP * h * k2.ds,
            e + HALF_STEP * h * k2.de,
            i + HALF_STEP * h * k2.di,
            r + HALF_STEP * h * k2.dr
        );

        // k4: slope at end of interval using k3
        const k4 = this.transitions.calculateDerivatives(
            s + h * k3.ds,
            e + h * k3.de,
            i + h * k3.di,
            r + h * k3.dr
        );

        // Helper to compute weighted average of slopes
        const weightedSum = (k1, k2, k3, k4) => 
            (h / WEIGHT_DIVISOR) * (k1 + WEIGHT_MIDDLE * k2 + WEIGHT_MIDDLE * k3 + k4);

        // Return new state using RK4 formula
        return [
            s + weightedSum(k1.ds, k2.ds, k3.ds, k4.ds),
            e + weightedSum(k1.de, k2.de, k3.de, k4.de),
            i + weightedSum(k1.di, k2.di, k3.di, k4.di),
            r + weightedSum(k1.dr, k2.dr, k3.dr, k4.dr)
        ];
    }
}

// ============================================================================
// SEIR Solver
// ============================================================================

/**
 * Main solver function for the SEIR model
 * 
 * @param {number} S0 - Initial proportion susceptible (0-1)
 * @param {number} R0 - Basic reproduction number
 * @param {number} latent_period - Days in exposed state
 * @param {number} infectious_period - Days in infectious state
 * @param {number} n - Number of time steps to simulate
 * @param {number} death_onset - Days until disease-induced death (default: 100)
 * @param {number} immunity_duration - Years of immunity (default: 1)
 * @param {number} life_expectancy - Years of life expectancy (default: 76)
 * @param {number} vaccination_rate - Proportion vaccinated at birth (0-1, default: 0.5)
 * @returns {Object} Time series data for S, E, I, R compartments
 */
export const solve = (
    S0,
    R0,
    latent_period,
    infectious_period,
    n,
    death_onset = 100,
    immunity_duration = 1,
    life_expectancy = 76,
    vaccination_rate = 0.5
) => {
    // Validate all input parameters
    ParameterValidator.validatePositive(n, 'Number of time steps');
    ParameterValidator.validateRange(S0, 0, 1, 'Initial susceptibility');
    ParameterValidator.validatePositive(R0, 'R0');
    ParameterValidator.validatePositive(infectious_period, 'Infectious period');
    ParameterValidator.validatePositive(latent_period, 'Latent period');
    ParameterValidator.validatePositive(immunity_duration, 'Immunity duration');
    ParameterValidator.validatePositive(life_expectancy, 'Life expectancy');
    ParameterValidator.validateNonNegative(death_onset, 'Death onset');
    ParameterValidator.validateRange(vaccination_rate, 0, 1, 'Vaccination rate');

    // Initialize state arrays (use Float64Array for performance)
    const state = {
        s: new Float64Array(n + 1),
        e: new Float64Array(n + 1),
        i: new Float64Array(n + 1),
        r: new Float64Array(n + 1)
    };

    // Set initial conditions
    state.s[0] = S0 - vaccination_rate;  // Susceptible minus those vaccinated
    state.e[0] = 1.0 - S0;               // Initially exposed
    state.i[0] = 0.0;                    // No infectious initially
    state.r[0] = vaccination_rate;       // Vaccinated individuals start as recovered

    // Setup model components
    const params = new SEIRParameters({
        R0, infectious_period, latent_period, death_onset, immunity_duration, life_expectancy
    });
    const transitions = new SEIRTransitions(params, vaccination_rate);
    const integrator = new RK4Integrator(transitions);

    // Integrate over time using RK4
    const timeStep = 1.0;  // 1 day per step
    for (let ix = 0; ix < n; ix++) {
        const [s_new, e_new, i_new, r_new] = integrator.step(
            [state.s[ix], state.e[ix], state.i[ix], state.r[ix]], 
            timeStep
        );

        // Clamp values to [0,1] to prevent numerical errors
        state.s[ix + 1] = clamp(s_new);
        state.e[ix + 1] = clamp(e_new);
        state.i[ix + 1] = clamp(i_new);
        state.r[ix + 1] = clamp(r_new);
    }

    // Return formatted data for plotting
    return {
        s: formatDataForPlot(state.s),
        e: formatDataForPlot(state.e),
        i: formatDataForPlot(state.i),
        r: formatDataForPlot(state.r),
        ymax: CONSTANTS.PERCENTAGE_SCALE
    };
};

// ============================================================================
// Plot Configuration
// ============================================================================

/**
 * Manages plot configuration and SVG setup
 */
class PlotConfiguration {
    constructor(plot_id, ctrl_id) {
        // Select DOM elements
        this.svg = d3.select(plot_id).append('svg');
        this.ctrls = d3.select(ctrl_id);
        
        // Plot layout configuration
        this.margin = { top: 20, right: 20, bottom: 30, left: 80 };
        this.axis_width = 2;
        
        // Default parameter values
        this.params = {
            S0: 0.99,                    // 99% initially susceptible
            R0: 3.0,                     // Basic reproduction number
            latent_period: 7.0,          // 7 days latent period
            infectious_period: 14.0,     // 14 days infectious period
            n_days: 3000,                // Simulate 3000 days
            y_max: 100,                  // Y-axis max (percentage)
            death_onset: 100,            // 100 days to death
            immunity_duration: 1,        // 1 year of immunity
            life_expectancy: 76,         // 76 years life expectancy
            vaccination_rate: 0.5        // 50% vaccination rate
        };

        this.updateDimensions();
        this.initializeLineGenerator();
    }

    /**
     * Updates plot dimensions from SVG bounding box
     */
    updateDimensions() {
        const svg_rect = this.svg.node().getBoundingClientRect();
        this.width = svg_rect.width;
        this.height = svg_rect.height;
    }

    /**
     * Initializes D3 line generator for drawing curves
     */
    initializeLineGenerator() {
        this.draw_line = d3.svg.line()
            .x(d => this.x_range(d.x))
            .y(d => this.y_range(d.y))
            .interpolate('linear');
    }
}

// ============================================================================
// Plot Parameter Management
// ============================================================================

/**
 * Manages parameter updates from form controls
 */
class PlotParameterManager {
    constructor(plot, param_vals) {
        this.plot = plot;
        this.param_vals = param_vals;  // Optional predefined parameter values
    }

    /**
     * Initializes parameters from form controls without updating plot
     */
    initialize() {
        this.plot.ctrls.selectAll('select').each(this.createSetter(false));
        this.plot.ctrls.selectAll('input').each(this.createSetter(false));
    }

    /**
     * Sets up event handlers for parameter changes
     */
    setupEventHandlers() {
        this.plot.ctrls.selectAll('select').on("change.param_val", this.createSetter(true));
        this.plot.ctrls.selectAll('input').on("change.param_val", this.createSetter(true));
        this.plot.ctrls.selectAll('input').on("input.param_val", this.createSetter(true));
    }

    /**
     * Creates a parameter setter function
     * @param {boolean} shouldUpdate - Whether to update plot after setting parameter
     * @returns {Function} Setter function to be called on form control
     */
    createSetter(shouldUpdate) {
        const plot = this.plot;
        const param_vals = this.param_vals;
        
        return function() {
            // Verify this control corresponds to a known parameter
            if (!(this.id in plot.params)) {
                console.log(`Form control for unknown parameter '${this.id}'`);
                return;
            }

            // Handle range inputs differently (they use discrete value arrays)
            if (this.type === "range") {
                PlotParameterManager.handleRangeInput(this, plot, param_vals);
            } else {
                plot.params[this.id] = parseFloat(this.value);
            }

            // Update plot if requested
            if (shouldUpdate) {
                plot.update();
            }
        };
    }

    /**
     * Handles range input controls with discrete value arrays
     */
    static handleRangeInput(input, plot, param_vals) {
        // Initialize range if not yet set up
        if (input.min === input.max) {
            if (!PlotParameterManager.initializeRangeInput(input, param_vals)) {
                return;
            }
        }

        // Get value from discrete array
        const ix = parseInt(input.value);
        plot.params[input.id] = input.data[ix].value;

        // Update displayed label
        const parent = d3.select(input.parentNode);
        const label = parent.select(".show_value")[0][0];
        label.textContent = input.data[ix].label !== undefined 
            ? input.data[ix].label 
            : input.data[ix].value;
    }

    /**
     * Initializes a range input with discrete values
     */
    static initializeRangeInput(input, param_vals) {
        if (!(input.id in param_vals)) {
            console.log(`No values to initialise '${input.id}'`);
            return false;
        }

        const values = param_vals[input.id];
        input.min = 0;
        input.max = values.length - 1;
        input.step = 1;
        input.data = values;

        // Set default value if specified
        const defaultIndex = values.findIndex(v => v.default);
        if (defaultIndex >= 0) {
            input.value = defaultIndex;
        }

        return true;
    }
}

// ============================================================================
// Plot Rendering
// ============================================================================

/**
 * Handles all plot rendering operations
 */
class PlotRenderer {
    constructor(plot) {
        this.plot = plot;
    }

    /**
     * Main update function - recalculates model and redraws plot
     */
    update() {
        // Solve SEIR model with current parameters
        const output = solve(
            this.plot.params.S0,
            this.plot.params.R0,
            this.plot.params.latent_period,
            this.plot.params.infectious_period,
            this.plot.params.n_days,
            this.plot.params.death_onset,
            this.plot.params.immunity_duration,
            this.plot.params.life_expectancy,
            this.plot.params.vaccination_rate
        );

        // Update plot components
        this.updateScales();
        this.updateAxes();
        this.drawDataSeries(output);
    }

    /**
     * Updates D3 scales and axes based on current dimensions and parameters
     */
    updateScales() {
        const plot = this.plot;

        // X-axis scale (time in days)
        if (!plot.x_range) {
            plot.x_range = d3.scale.linear();
        }
        plot.x_range
            .range([plot.margin.left, plot.width - plot.margin.right])
            .domain([0, plot.params.n_days]);

        // Y-axis scale (percentage)
        if (!plot.y_range) {
            plot.y_range = d3.scale.linear();
        }
        plot.y_range
            .range([plot.height - plot.margin.bottom, plot.margin.top])
            .domain([0, plot.params.y_max]);

        // X-axis generator
        if (!plot.x_axis) {
            plot.x_axis = d3.svg.axis();
        }
        plot.x_axis
            .scale(plot.x_range)
            .ticks(10)
            .tickSize(0);

        // Y-axis generator
        if (!plot.y_axis) {
            plot.y_axis = d3.svg.axis();
        }
        plot.y_axis
            .scale(plot.y_range)
            .orient('left')
            .tickSize(0)
            .tickFormat(d => `${d}%`)
            .ticks(4);
    }

    /**
     * Updates both X and Y axes
     */
    updateAxes() {
        this.updateXAxis();
        this.updateYAxis();
    }

    /**
     * Updates X-axis line and tick marks
     */
    updateXAxis() {
        const plot = this.plot;

        // X-axis baseline
        if (!plot.x_line) {
            plot.x_line = plot.svg.append('svg:line')
                .attr('class', 'x axis')
                .attr('stroke-width', plot.axis_width);
        }
        plot.x_line
            .attr("x1", plot.x_range(0) - plot.axis_width)
            .attr("y1", plot.y_range(0) + plot.axis_width)
            .attr("x2", plot.x_range(plot.params.n_days))
            .attr("y2", plot.y_range(0) + plot.axis_width);

        // X-axis tick marks and labels
        if (!plot.x_ticks) {
            plot.x_ticks = plot.svg.append('svg:g')
                .attr('class', 'x tick');
        }
        plot.x_ticks
            .attr('transform', `translate(0,${plot.height - 0.5 * plot.margin.bottom})`)
            .call(plot.x_axis);
    }

    /**
     * Updates Y-axis line and tick marks
     */
    updateYAxis() {
        const plot = this.plot;

        // Y-axis baseline
        if (!plot.y_line) {
            plot.y_line = plot.svg.append('svg:line')
                .attr('class', 'y axis')
                .attr('stroke-width', plot.axis_width);
        }
        plot.y_line
            .attr("x1", plot.x_range(0) - plot.axis_width)
            .attr("y1", plot.y_range(0) + plot.axis_width)
            .attr("x2", plot.x_range(0) - plot.axis_width)
            .attr("y2", plot.y_range(plot.params.y_max));

        // Y-axis tick marks and labels
        if (!plot.y_ticks) {
            plot.y_ticks = plot.svg.append('svg:g')
                .attr('class', 'y tick');
        }
        plot.y_ticks
            .attr('transform', `translate(${0.7 * plot.margin.left},0)`)
            .call(plot.y_axis);
    }

    /**
     * Draws the SEIR data series as SVG paths
     */
    drawDataSeries(output) {
        // Remove old series
        this.plot.svg.selectAll('path.series').remove();

        // Define series to draw (S, E, I, R)
        const series = [
            { data: output.s, class: 'varS' },  // Susceptible
            { data: output.e, class: 'varE' },  // Exposed
            { data: output.i, class: 'varI' },  // Infectious
            { data: output.r, class: 'varR' }   // Recovered
        ];

        // Draw each series as an SVG path
        series.forEach(({ data, class: className }) => {
            this.plot.svg.append('svg:path')
                .attr('d', this.plot.draw_line(data))
                .attr('class', `${className} series`);
        });
    }
}

// ============================================================================
// Main Plot Function
// ============================================================================

/**
 * Main entry point for creating an interactive SEIR plot
 * 
 * @param {string} plot_id - CSS selector for plot container
 * @param {string} ctrl_id - CSS selector for controls container
 * @param {Object} param_vals - Optional predefined parameter value arrays
 */
export const plot = (plot_id, ctrl_id, param_vals = {}) => {
    // Initialize plot components
    const plotConfig = new PlotConfiguration(plot_id, ctrl_id);
    const paramManager = new PlotParameterManager(plotConfig, param_vals);
    const renderer = new PlotRenderer(plotConfig);

    // Setup parameters and event handlers
    paramManager.initialize();
    paramManager.setupEventHandlers();

    // Attach update function to plot config
    plotConfig.update = () => renderer.update();

    // Handle window resize
    d3.select(window).on('resize', () => {
        plotConfig.updateDimensions();
        plotConfig.update();
    });

    // Initial render
    plotConfig.update();
};