// Enhanced SEIR model implementation with death rates, vaccination, and immunity waning
export const solve = (S0, R0, latent_period, infectious_period, n, 
                      death_onset = 100, immunity_duration = 1, life_expectancy = 76, 
                      vaccination_rate = 0.5) => {
    // Parameter validation
    if (n <= 0) throw new Error('Number of time steps must be positive');
    if (S0 < 0 || S0 > 1) throw new Error('Initial susceptibility must be between 0 and 1');
    if (R0 <= 0) throw new Error('R0 must be positive');
    if (infectious_period <= 0) throw new Error('Infectious period must be positive');
    if (latent_period <= 0) throw new Error('Latent period must be positive');
    if (immunity_duration <= 0) throw new Error('Immunity duration must be positive');
    if (life_expectancy <= 0) throw new Error('Life expectancy must be positive');
    if (death_onset < 0) throw new Error('Death onset must be positive'); // can be 0
    if (vaccination_rate < 0 || vaccination_rate > 1) throw new Error('Vaccination rate must be between 0 and 1');
    const s = new Float64Array(n + 1);
    const e = new Float64Array(n + 1);
    const i = new Float64Array(n + 1);
    const r = new Float64Array(n + 1);
    
    // Calculate derived parameters
    const alpha = death_onset > 0 ? 1 / death_onset : 0; // Death rate from infection
    const gamma = 1 / infectious_period; // Recovery rate
    const omega = 1 / (365 * immunity_duration); // Immunity waning rate
    const mu = 1 / (365 * life_expectancy); // Natural death rate
    const sigma = 1 / latent_period; // Latent to infectious rate
    const p = vaccination_rate; // Vaccination rate
    const beta = R0 * (gamma + mu + alpha) * (sigma + mu) / sigma; // Transmission rate
    
    // Initial conditions
    s[0] = S0 - p; // Account for vaccination
    e[0] = 1.0 - S0;
    i[0] = 0.0;
    r[0] = p; // Include vaccinated individuals

    // Enhanced ODE equations with death, vaccination, and immunity waning
    const s_to_e = (s, e, i, r) => beta * s * i;
    const e_to_i = (s, e, i, r) => sigma * e;
    const i_to_r = (s, e, i, r) => gamma * i;
    const r_to_s = (s, e, i, r) => omega * r;
    const s_death = (s, e, i, r) => mu * s;
    const e_death = (s, e, i, r) => mu * e;
    const i_death = (s, e, i, r) => (mu + alpha) * i;
    const r_death = (s, e, i, r) => mu * r;
    const birth_susceptible = (s, e, i, r) => mu * (1 - p);
    const birth_vaccinated = (s, e, i, r) => mu * p;

    // Runge-Kutta 4th order integration
    const rk4_step = (x, t, h) => {
        const [s_curr, e_curr, i_curr, r_curr] = x;
        
        // k1 = f(x, t)
        const k1_s = -s_to_e(s_curr, e_curr, i_curr, r_curr) + r_to_s(s_curr, e_curr, i_curr, r_curr) 
                    - s_death(s_curr, e_curr, i_curr, r_curr) + birth_susceptible(s_curr, e_curr, i_curr, r_curr);
        const k1_e = s_to_e(s_curr, e_curr, i_curr, r_curr) - e_to_i(s_curr, e_curr, i_curr, r_curr) 
                    - e_death(s_curr, e_curr, i_curr, r_curr);
        const k1_i = e_to_i(s_curr, e_curr, i_curr, r_curr) - i_to_r(s_curr, e_curr, i_curr, r_curr) 
                    - i_death(s_curr, e_curr, i_curr, r_curr);
        const k1_r = i_to_r(s_curr, e_curr, i_curr, r_curr) - r_to_s(s_curr, e_curr, i_curr, r_curr) 
                    - r_death(s_curr, e_curr, i_curr, r_curr) + birth_vaccinated(s_curr, e_curr, i_curr, r_curr);
        
        // k2 = f(x + 0.5*h*k1, t + 0.5*h)
        const x_k2 = [s_curr + 0.5*h*k1_s, e_curr + 0.5*h*k1_e, i_curr + 0.5*h*k1_i, r_curr + 0.5*h*k1_r];
        const k2_s = -s_to_e(...x_k2) + r_to_s(...x_k2) - s_death(...x_k2) + birth_susceptible(...x_k2);
        const k2_e = s_to_e(...x_k2) - e_to_i(...x_k2) - e_death(...x_k2);
        const k2_i = e_to_i(...x_k2) - i_to_r(...x_k2) - i_death(...x_k2);
        const k2_r = i_to_r(...x_k2) - r_to_s(...x_k2) - r_death(...x_k2) + birth_vaccinated(...x_k2);
        
        // k3 = f(x + 0.5*h*k2, t + 0.5*h)
        const x_k3 = [s_curr + 0.5*h*k2_s, e_curr + 0.5*h*k2_e, i_curr + 0.5*h*k2_i, r_curr + 0.5*h*k2_r];
        const k3_s = -s_to_e(...x_k3) + r_to_s(...x_k3) - s_death(...x_k3) + birth_susceptible(...x_k3);
        const k3_e = s_to_e(...x_k3) - e_to_i(...x_k3) - e_death(...x_k3);
        const k3_i = e_to_i(...x_k3) - i_to_r(...x_k3) - i_death(...x_k3);
        const k3_r = i_to_r(...x_k3) - r_to_s(...x_k3) - r_death(...x_k3) + birth_vaccinated(...x_k3);
        
        // k4 = f(x + h*k3, t + h)
        const x_k4 = [s_curr + h*k3_s, e_curr + h*k3_e, i_curr + h*k3_i, r_curr + h*k3_r];
        const k4_s = -s_to_e(...x_k4) + r_to_s(...x_k4) - s_death(...x_k4) + birth_susceptible(...x_k4);
        const k4_e = s_to_e(...x_k4) - e_to_i(...x_k4) - e_death(...x_k4);
        const k4_i = e_to_i(...x_k4) - i_to_r(...x_k4) - i_death(...x_k4);
        const k4_r = i_to_r(...x_k4) - r_to_s(...x_k4) - r_death(...x_k4) + birth_vaccinated(...x_k4);
        
        // Update using RK4 formula
        return [
            s_curr + (h/6) * (k1_s + 2*k2_s + 2*k3_s + k4_s),
            e_curr + (h/6) * (k1_e + 2*k2_e + 2*k3_e + k4_e),
            i_curr + (h/6) * (k1_i + 2*k2_i + 2*k3_i + k4_i),
            r_curr + (h/6) * (k1_r + 2*k2_r + 2*k3_r + k4_r)
        ];
    };

    // Integrate using RK4
    for (let ix = 0; ix < n; ix++) {
        const h = 1.0; // Time step
        const [s_new, e_new, i_new, r_new] = rk4_step([s[ix], e[ix], i[ix], r[ix]], ix, h);
        
        s[ix + 1] = Math.max(0, Math.min(1, s_new)); // Ensure values stay in [0,1]
        e[ix + 1] = Math.max(0, Math.min(1, e_new));
        i[ix + 1] = Math.max(0, Math.min(1, i_new));
        r[ix + 1] = Math.max(0, Math.min(1, r_new));
    }

    const data_S = [];
    const data_E = [];
    const data_I = [];
    const data_R = [];
    let data_max = 0.0;

    // Convert from population fractions to percentages.
    const scale_by = 100;

    for (let ix = 0; ix < n + 1; ix++) {
        data_S.push({x: ix, y: scale_by * s[ix]});
        data_E.push({x: ix, y: scale_by * e[ix]});
        data_I.push({x: ix, y: scale_by * i[ix]});
        data_R.push({x: ix, y: scale_by * r[ix]});
        
        if (data_S[ix].y > data_max) { data_max = data_S[ix].y; }
        if (data_E[ix].y > data_max) { data_max = data_E[ix].y; }
        if (data_I[ix].y > data_max) { data_max = data_I[ix].y; }
        if (data_R[ix].y > data_max) { data_max = data_R[ix].y; }
    }

    // Return lists of objects for use with D3.
    data_max = scale_by;
    return {s: data_S, e: data_E, i: data_I, r: data_R, ymax: data_max};
};

export const plot = (plot_id, ctrl_id, param_vals = {}) => {
    const plot = {};

    plot.svg = d3.select(plot_id).append('svg');

    const svg_rect = plot.svg.node().getBoundingClientRect();
    plot.width = svg_rect.width;
    plot.height = svg_rect.height;
    plot.margin = {
        top: 20,
        right: 20,
        bottom: 30,
        left: 80
    };
    plot.axis_width = 2;

    plot.ctrls = d3.select(ctrl_id);

    plot.params = {
        S0: 0.99,
        R0: 3.0,
        latent_period: 7.0, // days
        infectious_period: 14.0, // days
        n_days: 3000, // Extended simulation period
        y_max: 100,
        death_onset: 100, // days
        immunity_duration: 1, // years
        life_expectancy: 76, // years
        vaccination_rate: 0.5 // 50% vaccination rate
    };

    // Set parameters to initial form values.
    const set_param = (update_plot) => {
        return function() {
            if (this.id in plot.params) {
                if (this.type === "range") {
                    if (this.min === this.max) {
                        if (this.id in param_vals) {
                            const values = param_vals[this.id];
                            this.min = 0;
                            this.max = values.length - 1;
                            this.step = 1;
                            this.data = values;
                            // Pick the default initial value, if specified.
                            const def = values.findIndex(v => v.default);
                            if (def >= 0) {
                                this.value = def;
                            }
                        } else {
                            console.log(`No values to initialise '${this.id}'`);
                            return;
                        }
                    }
                    const ix = parseInt(this.value);
                    // Update the parameter value.
                    plot.params[this.id] = this.data[ix].value;
                    // Display the parameter value.
                    const parent = d3.select(this.parentNode);
                    const label = parent.select(".show_value")[0][0];
                    if (this.data[ix].label !== undefined) {
                        label.textContent = this.data[ix].label;
                    } else {
                        label.textContent = this.data[ix].value;
                    }
                } else {
                    plot.params[this.id] = parseFloat(this.value);
                }
                if (update_plot) {
                    plot.update();
                }
            } else {
                console.log(`Form control for unknown parameter '${this.id}'`);
            }
        };
    };

    plot.ctrls.selectAll('select').each(set_param(false));
    plot.ctrls.selectAll('input').each(set_param(false));

    plot.draw_line = d3.svg.line()
        .x(d => plot.x_range(d.x))
        .y(d => plot.y_range(d.y))
        .interpolate('linear');

    plot.update = () => {
        const output = solve(
            plot.params.S0,
            plot.params.R0,
            plot.params.latent_period,
            plot.params.infectious_period,
            plot.params.n_days,
            plot.params.death_onset,
            plot.params.immunity_duration,
            plot.params.life_expectancy,
            plot.params.vaccination_rate);

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

        // (Re)draw axes.
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
        if (plot.x_ticks === undefined) {
            plot.x_ticks = plot.svg.append('svg:g')
                .attr('class', 'x tick');
        }
        plot.x_ticks
            .attr('transform', `translate(0,${plot.height - 0.5 * plot.margin.bottom})`)
            .call(plot.x_axis);
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
        if (plot.y_ticks === undefined) {
            plot.y_ticks = plot.svg.append('svg:g')
                .attr('class', 'y tick');
        }
        plot.y_ticks
            .attr('transform', `translate(${0.7 * plot.margin.left},0)`)
            .call(plot.y_axis);

        // Remove existing data series.
        plot.svg.selectAll('path.series').remove();

        plot.svg.append('svg:path')
            .attr('d', plot.draw_line(output.s))
            .attr('class', 'varS series');
        plot.svg.append('svg:path')
            .attr('d', plot.draw_line(output.e))
            .attr('class', 'varE series');
        plot.svg.append('svg:path')
            .attr('d', plot.draw_line(output.i))
            .attr('class', 'varI series');
        plot.svg.append('svg:path')
            .attr('d', plot.draw_line(output.r))
            .attr('class', 'varR series');
    };

    // Add update handlers for each input element.
    plot.ctrls.selectAll('select').on("change.param_val", set_param(true));
    plot.ctrls.selectAll('input').on("change.param_val", set_param(true));
    plot.ctrls.selectAll('input').on("input.param_val", set_param(true));

    d3.select(window).on('resize', () => {
        const svg_rect = plot.svg.node().getBoundingClientRect();
        plot.width = svg_rect.width;
        plot.height = svg_rect.height;
        plot.update();
    });

    plot.update();
};
