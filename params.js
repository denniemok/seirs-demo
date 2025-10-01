// Parameter values for the SEIRS model
// Helper function to generate parameter arrays
const generateParams = (start, end, step, defaultVal = null, labelFunc = null) => {
    const params = [];
    // Use integer arithmetic to avoid floating-point precision issues
    const numSteps = Math.round((end - start) / step);
    
    for (let i = 0; i <= numSteps; i++) {
        const value = Math.round((start + i * step) * 100) / 100;
        const param = { value };
        
        if (value === defaultVal) param.default = true;
        if (labelFunc) param.label = labelFunc(value);
        
        params.push(param);
    }
    return params;
};

export const param_vals = {
    n_days: generateParams(50, 3000, 50, 3000),
    y_max: generateParams(5, 100, 5, 100),
    R0: generateParams(1, 5, 0.1, 3.0),
    latent_period: generateParams(1, 30, 1, 7, val => `${val} day${val !== 1 ? 's' : ''}`),
    infectious_period: generateParams(1, 30, 1, 14, val => `${val} day${val !== 1 ? 's' : ''}`),
    S0: generateParams(0, 1, 0.01, 0.99, val => `${Math.round(val * 100)}%`),
    death_onset: generateParams(0, 1000, 1, 15, val => `${val} day${val !== 1 ? 's' : ''}`),
    immunity_duration: generateParams(0.1, 10, 0.1, 1.0, val => `${val} year${val !== 1 ? 's' : ''}`),
    life_expectancy: generateParams(1, 100, 1, 76, val => `${val} year${val !== 1 ? 's' : ''}`),
    vaccination_rate: generateParams(0, 1, 0.01, 0.0, val => `${Math.round(val * 100)}%`)
};
