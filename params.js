// Parameter values for the SEIRS model
// Helper function to generate parameter arrays
const generateParams = (start, end, step, defaultVal = null, labelFunc = null) => {
    const params = [];
    // Use integer arithmetic to avoid floating-point precision issues
    const numSteps = Math.round((end - start) / step);
    
    for (let i = 0; i <= numSteps; i++) {
        const value = Math.round((start + i * step) * 100) / 100; // Round to 2 decimal places
        const param = { value };
        
        if (value === defaultVal) {
            param.default = true;
        }
        
        if (labelFunc) {
            param.label = labelFunc(value);
        }
        
        params.push(param);
    }
    return params;
};

export const param_vals = {
    // n_days: from 50 to 3000, interval 50
    n_days: generateParams(50, 3000, 50, 3000),
    
    // y_max: 5 to 100, interval 5
    y_max: generateParams(5, 100, 5, 100),
    
    // R0: 1 to 5, interval 0.1
    R0: generateParams(1, 5, 0.1, 3.0),
    
    // latent_period: 1 to 30, interval 1
    latent_period: generateParams(1, 30, 1, 7, (val) => `${val} day${val !== 1 ? 's' : ''}`),
    
    // infectious_period: 1 to 30, interval 1
    infectious_period: generateParams(1, 30, 1, 14, (val) => `${val} day${val !== 1 ? 's' : ''}`),
    
    // S0: from 0 to 1, interval 0.01
    S0: generateParams(0, 1, 0.01, 0.99, (val) => `${Math.round(val * 100)}%`),
    
    // death_onset: 0 to 1000, interval 1
    death_onset: generateParams(0, 1000, 1, 15, (val) => `${val} day${val !== 1 ? 's' : ''}`),
    
    // immunity_duration: 0.1 to 10, interval 0.1
    immunity_duration: generateParams(0.1, 10, 0.1, 1.0, (val) => `${val} year${val !== 1 ? 's' : ''}`),
    
    // life_expectancy: 1 to 100, interval 1
    life_expectancy: generateParams(1, 100, 1, 76, (val) => `${val} year${val !== 1 ? 's' : ''}`),
    
    // vaccination_rate: 0 to 1, interval 0.01
    vaccination_rate: generateParams(0, 1, 0.01, 0.0, (val) => `${Math.round(val * 100)}%`)
};
