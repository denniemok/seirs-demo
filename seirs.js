// Enhanced SEIR model implementation with death rates, vaccination, and immunity waning

// Constants
const DAYS_PER_YEAR = 365;
const PERCENTAGE_SCALE = 100;
const RK4_WEIGHT_MIDDLE = 2;
const RK4_WEIGHT_DIVISOR = 6;
const RK4_HALF_STEP = 0.5;

// Validation helpers
const validateParameter = (value, min, max, name) => {
    if (value < min || value > max) {
        throw new Error(`${name} must be between ${min} and ${max}`);
    }
};

const validatePositive = (value, name) => {
    if (value <= 0) {
        throw new Error(`${name} must be positive`);
    }
};

const validateNonNegative = (value, name) => {
    if (value < 0) {
        throw new Error(`${name} must be non-negative`);
    }
};

// Model parameters calculator
class SEIRParameters {
    constructor(R0, infectious_period, latent_period, death_onset, immunity_duration, life_expectancy) {
        this.alpha = death_onset > 0 ? 1 / death_onset : 0; // Death rate from infection
        this.gamma = 1 / infectious_period; // Recovery rate
        this.omega = 1 / (DAYS_PER_YEAR * immunity_duration); // Immunity waning rate
        this.mu = 1 / (DAYS_PER_YEAR * life_expectancy); // Natural death rate
        this.sigma = 1 / latent_period; // Latent to infectious rate
        this.beta = R0 * (this.gamma + this.mu + this.alpha) * (this.sigma + this.mu) / this.sigma; // Transmission rate
    }
}

// SEIR state transitions
class SEIRTransitions {
    constructor(params, vaccination_rate) {
        this.params = params;
        this.p = vaccination_rate;
    }

    susceptibleToExposed(s, i) {
        return this.params.beta * s * i;
    }

    exposedToInfectious(e) {
        return this.params.sigma * e;
    }

    infectiousToRecovered(i) {
        return this.params.gamma * i;
    }

    recoveredToSusceptible(r) {
        return this.params.omega * r;
    }

    susceptibleDeath(s) {
        return this.params.mu * s;
    }

    exposedDeath(e) {
        return this.params.mu * e;
    }

    infectiousDeath(i) {
        return (this.params.mu + this.params.alpha) * i;
    }

    recoveredDeath(r) {
        return this.params.mu * r;
    }

    birthSusceptible() {
        return this.params.mu * (1 - this.p);
    }

    birthVaccinated() {
        return this.params.mu * this.p;
    }

    // Calculate derivatives for each compartment
    calculateDerivatives(s, e, i, r) {
        const s_to_e = this.susceptibleToExposed(s, i);
        const e_to_i = this.exposedToInfectious(e);
        const i_to_r = this.infectiousToRecovered(i);
        const r_to_s = this.recoveredToSusceptible(r);

        return {
            ds: -s_to_e + r_to_s - this.susceptibleDeath(s) + this.birthSusceptible(),
            de: s_to_e - e_to_i - this.exposedDeath(e),
            di: e_to_i - i_to_r - this.infectiousDeath(i),
            dr: i_to_r - r_to_s - this.recoveredDeath(r) + this.birthVaccinated()
        };
    }
}

// Runge-Kutta 4th order integrator
class RK4Integrator {
    constructor(transitions) {
        this.transitions = transitions;
    }

    step(state, h) {
        const [s, e, i, r] = state;

        // k1 = f(x, t)
        const k1 = this.transitions.calculateDerivatives(s, e, i, r);

        // k2 = f(x + 0.5*h*k1, t + 0.5*h)
        const k2 = this.transitions.calculateDerivatives(
            s + RK4_HALF_STEP * h * k1.ds,
            e + RK4_HALF_STEP * h * k1.de,
            i + RK4_HALF_STEP * h * k1.di,
            r + RK4_HALF_STEP * h * k1.dr
        );

        // k3 = f(x + 0.5*h*k2, t + 0.5*h)
        const k3 = this.transitions.calculateDerivatives(
            s + RK4_HALF_STEP * h * k2.ds,
            e + RK4_HALF_STEP * h * k2.de,
            i + RK4_HALF_STEP * h * k2.di,
            r + RK4_HALF_STEP * h * k2.dr
        );

        // k4 = f(x + h*k3, t + h)
        const k4 = this.transitions.calculateDerivatives(
            s + h * k3.ds,
            e + h * k3.de,
            i + h * k3.di,
            r + h * k3.dr
        );

        // Update using RK4 formula: x_new = x + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        return [
            s + (h / RK4_WEIGHT_DIVISOR) * (k1.ds + RK4_WEIGHT_MIDDLE * k2.ds + RK4_WEIGHT_MIDDLE * k3.ds + k4.ds),
            e + (h / RK4_WEIGHT_DIVISOR) * (k1.de + RK4_WEIGHT_MIDDLE * k2.de + RK4_WEIGHT_MIDDLE * k3.de + k4.de),
            i + (h / RK4_WEIGHT_DIVISOR) * (k1.di + RK4_WEIGHT_MIDDLE * k2.di + RK4_WEIGHT_MIDDLE * k3.di + k4.di),
            r + (h / RK4_WEIGHT_DIVISOR) * (k1.dr + RK4_WEIGHT_MIDDLE * k2.dr + RK4_WEIGHT_MIDDLE * k3.dr + k4.dr)
        ];
    }
}

// Clamp value between 0 and 1
const clamp = (value) => Math.max(0, Math.min(1, value));

// Main solver function
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
    // Parameter validation
    validatePositive(n, 'Number of time steps');
    validateParameter(S0, 0, 1, 'Initial susceptibility');
    validatePositive(R0, 'R0');
    validatePositive(infectious_period, 'Infectious period');
    validatePositive(latent_period, 'Latent period');
    validatePositive(immunity_duration, 'Immunity duration');
    validatePositive(life_expectancy, 'Life expectancy');
    validateNonNegative(death_onset, 'Death onset');
    validateParameter(vaccination_rate, 0, 1, 'Vaccination rate');

    // Initialize state arrays
    const s = new Float64Array(n + 1);
    const e = new Float64Array(n + 1);
    const i = new Float64Array(n + 1);
    const r = new Float64Array(n + 1);

    // Calculate model parameters
    const params = new SEIRParameters(R0, infectious_period, latent_period, death_onset, immunity_duration, life_expectancy);
    const transitions = new SEIRTransitions(params, vaccination_rate);
    const integrator = new RK4Integrator(transitions);

    // Set initial conditions
    s[0] = S0 - vaccination_rate; // Account for vaccination
    e[0] = 1.0 - S0;
    i[0] = 0.0;
    r[0] = vaccination_rate; // Include vaccinated individuals

    // Integrate using RK4
    const timeStep = 1.0;
    for (let ix = 0; ix < n; ix++) {
        const [s_new, e_new, i_new, r_new] = integrator.step([s[ix], e[ix], i[ix], r[ix]], timeStep);

        // Ensure values stay in [0,1]
        s[ix + 1] = clamp(s_new);
        e[ix + 1] = clamp(e_new);
        i[ix + 1] = clamp(i_new);
        r[ix + 1] = clamp(r_new);
    }

    // Convert results to output format
    const formatData = (array) => {
        return Array.from(array, (value, index) => ({
            x: index,
            y: PERCENTAGE_SCALE * value
        }));
    };

    return {
        s: formatData(s),
        e: formatData(e),
        i: formatData(i),
        r: formatData(r),
        ymax: PERCENTAGE_SCALE
    };
};

// Plot configuration and management
export const plot = (plot_id, ctrl_id, param_vals = {}) => {
    const plot = createPlotObject(plot_id, ctrl_id);
    
    initializePlotParameters(plot, param_vals);
    setupEventHandlers(plot);
    
    plot.update();
};

// Create plot object with initial configuration
const createPlotObject = (plot_id, ctrl_id) => {
    const plot = {
        svg: d3.select(plot_id).append('svg'),
        ctrls: d3.select(ctrl_id),
        margin: { top: 20, right: 20, bottom: 30, left: 80 },
        axis_width: 2,
        params: {
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
        }
    };

    const svg_rect = plot.svg.node().getBoundingClientRect();
    plot.width = svg_rect.width;
    plot.height = svg_rect.height;

    plot.draw_line = d3.svg.line()
        .x(d => plot.x_range(d.x))
        .y(d => plot.y_range(d.y))
        .interpolate('linear');

    plot.update = () => updatePlot(plot);

    return plot;
};

// Initialize plot parameters from form controls
const initializePlotParameters = (plot, param_vals) => {
    const setParam = createParamSetter(plot, param_vals);
    
    plot.ctrls.selectAll('select').each(setParam(false));
    plot.ctrls.selectAll('input').each(setParam(false));
};

// Create parameter setter function
const createParamSetter = (plot, param_vals) => {
    return (update_plot) => {
        return function() {
            if (!(this.id in plot.params)) {
                console.log(`Form control for unknown parameter '${this.id}'`);
                return;
            }

            if (this.type === "range") {
                handleRangeInput(this, plot, param_vals);
            } else {
                plot.params[this.id] = parseFloat(this.value);
            }

            if (update_plot) {
                plot.update();
            }
        };
    };
};

// Handle range input controls
const handleRangeInput = (input, plot, param_vals) => {
    if (input.min === input.max) {
        if (!initializeRangeInput(input, param_vals)) {
            return;
        }
    }

    const ix = parseInt(input.value);
    plot.params[input.id] = input.data[ix].value;

    // Display the parameter value
    const parent = d3.select(input.parentNode);
    const label = parent.select(".show_value")[0][0];
    label.textContent = input.data[ix].label !== undefined 
        ? input.data[ix].label 
        : input.data[ix].value;
};

// Initialize range input with values
const initializeRangeInput = (input, param_vals) => {
    if (!(input.id in param_vals)) {
        console.log(`No values to initialise '${input.id}'`);
        return false;
    }

    const values = param_vals[input.id];
    input.min = 0;
    input.max = values.length - 1;
    input.step = 1;
    input.data = values;

    // Pick the default initial value, if specified
    const def = values.findIndex(v => v.default);
    if (def >= 0) {
        input.value = def;
    }

    return true;
};

// Setup event handlers for plot controls
const setupEventHandlers = (plot) => {
    const setParam = createParamSetter(plot, {});
    
    plot.ctrls.selectAll('select').on("change.param_val", setParam(true));
    plot.ctrls.selectAll('input').on("change.param_val", setParam(true));
    plot.ctrls.selectAll('input').on("input.param_val", setParam(true));

    d3.select(window).on('resize', () => {
        const svg_rect = plot.svg.node().getBoundingClientRect();
        plot.width = svg_rect.width;
        plot.height = svg_rect.height;
        plot.update();
    });
};

// Update plot with current parameters
const updatePlot = (plot) => {
    const output = solve(
        plot.params.S0,
        plot.params.R0,
        plot.params.latent_period,
        plot.params.infectious_period,
        plot.params.n_days,
        plot.params.death_onset,
        plot.params.immunity_duration,
        plot.params.life_expectancy,
        plot.params.vaccination_rate
    );

    updateScales(plot);
    updateAxes(plot);
    drawDataSeries(plot, output);
};

// Update plot scales
const updateScales = (plot) => {
    if (plot.x_range === undefined) {
        plot.x_range = d3.scale.linear();
    }
    plot.x_range
        .range([plot.margin.left, plot.width - plot.margin.right])
        .domain([0, plot.params.n_days]);

    if (plot.y_range === undefined) {
        plot.y_range = d3.scale.linear();
    }
    plot.y_range
        .range([plot.height - plot.margin.bottom, plot.margin.top])
        .domain([0, plot.params.y_max]);

    if (plot.x_axis === undefined) {
        plot.x_axis = d3.svg.axis();
    }
    plot.x_axis
        .scale(plot.x_range)
        .ticks(10)
        .tickSize(0);

    if (plot.y_axis === undefined) {
        plot.y_axis = d3.svg.axis();
    }
    plot.y_axis
        .scale(plot.y_range)
        .orient('left')
        .tickSize(0)
        .tickFormat(d => `${d}%`)
        .ticks(4);
};

// Update plot axes
const updateAxes = (plot) => {
    // X-axis line
    if (plot.x_line === undefined) {
        plot.x_line = plot.svg.append('svg:line')
            .attr('class', 'x axis')
            .attr('stroke-width', plot.axis_width);
    }
    plot.x_line
        .attr("x1", plot.x_range(0) - plot.axis_width)
        .attr("y1", plot.y_range(0) + plot.axis_width)
        .attr("x2", plot.x_range(plot.params.n_days))
        .attr("y2", plot.y_range(0) + plot.axis_width);

    // X-axis ticks
    if (plot.x_ticks === undefined) {
        plot.x_ticks = plot.svg.append('svg:g')
            .attr('class', 'x tick');
    }
    plot.x_ticks
        .attr('transform', `translate(0,${plot.height - 0.5 * plot.margin.bottom})`)
        .call(plot.x_axis);

    // Y-axis line
    if (plot.y_line === undefined) {
        plot.y_line = plot.svg.append('svg:line')
            .attr('class', 'y axis')
            .attr('stroke-width', plot.axis_width);
    }
    plot.y_line
        .attr("x1", plot.x_range(0) - plot.axis_width)
        .attr("y1", plot.y_range(0) + plot.axis_width)
        .attr("x2", plot.x_range(0) - plot.axis_width)
        .attr("y2", plot.y_range(plot.params.y_max));

    // Y-axis ticks
    if (plot.y_ticks === undefined) {
        plot.y_ticks = plot.svg.append('svg:g')
            .attr('class', 'y tick');
    }
    plot.y_ticks
        .attr('transform', `translate(${0.7 * plot.margin.left},0)`)
        .call(plot.y_axis);
};

// Draw data series on plot
const drawDataSeries = (plot, output) => {
    // Remove existing data series
    plot.svg.selectAll('path.series').remove();

    // Draw new data series
    const series = [
        { data: output.s, class: 'varS' },
        { data: output.e, class: 'varE' },
        { data: output.i, class: 'varI' },
        { data: output.r, class: 'varR' }
    ];

    series.forEach(({ data, class: className }) => {
        plot.svg.append('svg:path')
            .attr('d', plot.draw_line(data))
            .attr('class', `${className} series`);
    });
};