// Enhanced SEIR model implementation with death rates, vaccination, and immunity waning

// ============================================================================
// Constants
// ============================================================================

const CONSTANTS = {
    DAYS_PER_YEAR: 365,
    PERCENTAGE_SCALE: 100,
    RK4: {
        WEIGHT_MIDDLE: 2,
        WEIGHT_DIVISOR: 6,
        HALF_STEP: 0.5
    }
};

// ============================================================================
// Validation Utilities
// ============================================================================

class ParameterValidator {
    static validateRange(value, min, max, name) {
        if (value < min || value > max) {
            throw new Error(`${name} must be between ${min} and ${max}`);
        }
    }

    static validatePositive(value, name) {
        if (value <= 0) {
            throw new Error(`${name} must be positive`);
        }
    }

    static validateNonNegative(value, name) {
        if (value < 0) {
            throw new Error(`${name} must be non-negative`);
        }
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

const clamp = (value, min = 0, max = 1) => Math.max(min, Math.min(max, value));

const formatDataForPlot = (array, scale = CONSTANTS.PERCENTAGE_SCALE) => {
    return Array.from(array, (value, index) => ({
        x: index,
        y: scale * value
    }));
};

// ============================================================================
// SEIR Model Core
// ============================================================================

class SEIRParameters {
    constructor(config) {
        const { R0, infectious_period, latent_period, death_onset, immunity_duration, life_expectancy } = config;
        
        this.alpha = death_onset > 0 ? 1 / death_onset : 0;
        this.gamma = 1 / infectious_period;
        this.omega = 1 / (CONSTANTS.DAYS_PER_YEAR * immunity_duration);
        this.mu = 1 / (CONSTANTS.DAYS_PER_YEAR * life_expectancy);
        this.sigma = 1 / latent_period;
        this.beta = R0 * (this.gamma + this.mu + this.alpha) * (this.sigma + this.mu) / this.sigma;
    }
}

class SEIRTransitions {
    constructor(params, vaccination_rate) {
        this.params = params;
        this.vaccination_rate = vaccination_rate;
    }

    calculateDerivatives(s, e, i, r) {
        const { beta, sigma, gamma, omega, mu, alpha } = this.params;
        const p = this.vaccination_rate;

        const s_to_e = beta * s * i;
        const e_to_i = sigma * e;
        const i_to_r = gamma * i;
        const r_to_s = omega * r;

        return {
            ds: -s_to_e + r_to_s - mu * s + mu * (1 - p),
            de: s_to_e - e_to_i - mu * e,
            di: e_to_i - i_to_r - (mu + alpha) * i,
            dr: i_to_r - r_to_s - mu * r + mu * p
        };
    }
}

class RK4Integrator {
    constructor(transitions) {
        this.transitions = transitions;
    }

    step(state, h) {
        const [s, e, i, r] = state;
        const { HALF_STEP, WEIGHT_MIDDLE, WEIGHT_DIVISOR } = CONSTANTS.RK4;

        const k1 = this.transitions.calculateDerivatives(s, e, i, r);
        
        const k2 = this.transitions.calculateDerivatives(
            s + HALF_STEP * h * k1.ds,
            e + HALF_STEP * h * k1.de,
            i + HALF_STEP * h * k1.di,
            r + HALF_STEP * h * k1.dr
        );

        const k3 = this.transitions.calculateDerivatives(
            s + HALF_STEP * h * k2.ds,
            e + HALF_STEP * h * k2.de,
            i + HALF_STEP * h * k2.di,
            r + HALF_STEP * h * k2.dr
        );

        const k4 = this.transitions.calculateDerivatives(
            s + h * k3.ds,
            e + h * k3.de,
            i + h * k3.di,
            r + h * k3.dr
        );

        const weightedSum = (k1, k2, k3, k4) => 
            (h / WEIGHT_DIVISOR) * (k1 + WEIGHT_MIDDLE * k2 + WEIGHT_MIDDLE * k3 + k4);

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
    // Validate all parameters
    ParameterValidator.validatePositive(n, 'Number of time steps');
    ParameterValidator.validateRange(S0, 0, 1, 'Initial susceptibility');
    ParameterValidator.validatePositive(R0, 'R0');
    ParameterValidator.validatePositive(infectious_period, 'Infectious period');
    ParameterValidator.validatePositive(latent_period, 'Latent period');
    ParameterValidator.validatePositive(immunity_duration, 'Immunity duration');
    ParameterValidator.validatePositive(life_expectancy, 'Life expectancy');
    ParameterValidator.validateNonNegative(death_onset, 'Death onset');
    ParameterValidator.validateRange(vaccination_rate, 0, 1, 'Vaccination rate');

    // Initialize state arrays
    const state = {
        s: new Float64Array(n + 1),
        e: new Float64Array(n + 1),
        i: new Float64Array(n + 1),
        r: new Float64Array(n + 1)
    };

    // Set initial conditions
    state.s[0] = S0 - vaccination_rate;
    state.e[0] = 1.0 - S0;
    state.i[0] = 0.0;
    state.r[0] = vaccination_rate;

    // Setup model components
    const params = new SEIRParameters({
        R0, infectious_period, latent_period, death_onset, immunity_duration, life_expectancy
    });
    const transitions = new SEIRTransitions(params, vaccination_rate);
    const integrator = new RK4Integrator(transitions);

    // Integrate over time
    const timeStep = 1.0;
    for (let ix = 0; ix < n; ix++) {
        const [s_new, e_new, i_new, r_new] = integrator.step(
            [state.s[ix], state.e[ix], state.i[ix], state.r[ix]], 
            timeStep
        );

        state.s[ix + 1] = clamp(s_new);
        state.e[ix + 1] = clamp(e_new);
        state.i[ix + 1] = clamp(i_new);
        state.r[ix + 1] = clamp(r_new);
    }

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

class PlotConfiguration {
    constructor(plot_id, ctrl_id) {
        this.svg = d3.select(plot_id).append('svg');
        this.ctrls = d3.select(ctrl_id);
        this.margin = { top: 20, right: 20, bottom: 30, left: 80 };
        this.axis_width = 2;
        
        this.params = {
            S0: 0.99,
            R0: 3.0,
            latent_period: 7.0,
            infectious_period: 14.0,
            n_days: 3000,
            y_max: 100,
            death_onset: 100,
            immunity_duration: 1,
            life_expectancy: 76,
            vaccination_rate: 0.5
        };

        this.updateDimensions();
        this.initializeLineGenerator();
    }

    updateDimensions() {
        const svg_rect = this.svg.node().getBoundingClientRect();
        this.width = svg_rect.width;
        this.height = svg_rect.height;
    }

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

class PlotParameterManager {
    constructor(plot, param_vals) {
        this.plot = plot;
        this.param_vals = param_vals;
    }

    initialize() {
        this.plot.ctrls.selectAll('select').each(this.createSetter(false));
        this.plot.ctrls.selectAll('input').each(this.createSetter(false));
    }

    setupEventHandlers() {
        this.plot.ctrls.selectAll('select').on("change.param_val", this.createSetter(true));
        this.plot.ctrls.selectAll('input').on("change.param_val", this.createSetter(true));
        this.plot.ctrls.selectAll('input').on("input.param_val", this.createSetter(true));
    }

    createSetter(shouldUpdate) {
        const plot = this.plot;
        const param_vals = this.param_vals;
        
        return function() {
            if (!(this.id in plot.params)) {
                console.log(`Form control for unknown parameter '${this.id}'`);
                return;
            }

            if (this.type === "range") {
                PlotParameterManager.handleRangeInput(this, plot, param_vals);
            } else {
                plot.params[this.id] = parseFloat(this.value);
            }

            if (shouldUpdate) {
                plot.update();
            }
        };
    }

    static handleRangeInput(input, plot, param_vals) {
        if (input.min === input.max) {
            if (!PlotParameterManager.initializeRangeInput(input, param_vals)) {
                return;
            }
        }

        const ix = parseInt(input.value);
        plot.params[input.id] = input.data[ix].value;

        const parent = d3.select(input.parentNode);
        const label = parent.select(".show_value")[0][0];
        label.textContent = input.data[ix].label !== undefined 
            ? input.data[ix].label 
            : input.data[ix].value;
    }

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

class PlotRenderer {
    constructor(plot) {
        this.plot = plot;
    }

    update() {
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

        this.updateScales();
        this.updateAxes();
        this.drawDataSeries(output);
    }

    updateScales() {
        const plot = this.plot;

        if (!plot.x_range) {
            plot.x_range = d3.scale.linear();
        }
        plot.x_range
            .range([plot.margin.left, plot.width - plot.margin.right])
            .domain([0, plot.params.n_days]);

        if (!plot.y_range) {
            plot.y_range = d3.scale.linear();
        }
        plot.y_range
            .range([plot.height - plot.margin.bottom, plot.margin.top])
            .domain([0, plot.params.y_max]);

        if (!plot.x_axis) {
            plot.x_axis = d3.svg.axis();
        }
        plot.x_axis
            .scale(plot.x_range)
            .ticks(10)
            .tickSize(0);

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

    updateAxes() {
        this.updateXAxis();
        this.updateYAxis();
    }

    updateXAxis() {
        const plot = this.plot;

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

        if (!plot.x_ticks) {
            plot.x_ticks = plot.svg.append('svg:g')
                .attr('class', 'x tick');
        }
        plot.x_ticks
            .attr('transform', `translate(0,${plot.height - 0.5 * plot.margin.bottom})`)
            .call(plot.x_axis);
    }

    updateYAxis() {
        const plot = this.plot;

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

        if (!plot.y_ticks) {
            plot.y_ticks = plot.svg.append('svg:g')
                .attr('class', 'y tick');
        }
        plot.y_ticks
            .attr('transform', `translate(${0.7 * plot.margin.left},0)`)
            .call(plot.y_axis);
    }

    drawDataSeries(output) {
        this.plot.svg.selectAll('path.series').remove();

        const series = [
            { data: output.s, class: 'varS' },
            { data: output.e, class: 'varE' },
            { data: output.i, class: 'varI' },
            { data: output.r, class: 'varR' }
        ];

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

export const plot = (plot_id, ctrl_id, param_vals = {}) => {
    const plotConfig = new PlotConfiguration(plot_id, ctrl_id);
    const paramManager = new PlotParameterManager(plotConfig, param_vals);
    const renderer = new PlotRenderer(plotConfig);

    paramManager.initialize();
    paramManager.setupEventHandlers();

    plotConfig.update = () => renderer.update();

    d3.select(window).on('resize', () => {
        plotConfig.updateDimensions();
        plotConfig.update();
    });

    plotConfig.update();
};