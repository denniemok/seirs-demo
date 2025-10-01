// Parameter values for the SEIRS model

/**
 * Generates parameter configurations for model inputs
 * @param {number} start - Starting value
 * @param {number} end - Ending value
 * @param {number} step - Increment step
 * @param {number|null} defaultVal - Default value to mark
 * @param {Function|null} labelFunc - Optional label formatter
 * @returns {Array} Array of parameter objects
 */
const generateParams = (start, end, step, defaultVal = null, labelFunc = null) => {
    const numSteps = Math.round((end - start) / step);
    const params = [];
    
    for (let i = 0; i <= numSteps; i++) {
        const value = roundToTwoDecimals(start + i * step);
        const param = createParam(value, defaultVal, labelFunc);
        params.push(param);
    }
    
    return params;
};

const roundToTwoDecimals = (num) => Math.round(num * 100) / 100;

const createParam = (value, defaultVal, labelFunc) => {
    const param = { value };
    
    if (value === defaultVal) {
        param.default = true;
    }
    
    if (labelFunc) {
        param.label = labelFunc(value);
    }
    
    return param;
};

// Label formatters
const dayLabel = (val) => `${val} day${val !== 1 ? 's' : ''}`;
const yearLabel = (val) => `${val} year${val !== 1 ? 's' : ''}`;
const percentLabel = (val) => `${Math.round(val * 100)}%`;

export const param_vals = {
    n_days: generateParams(50, 3000, 50, 3000),
    y_max: generateParams(5, 100, 5, 100),
    R0: generateParams(1, 5, 0.1, 3.0),
    latent_period: generateParams(1, 30, 1, 7, dayLabel),
    infectious_period: generateParams(1, 30, 1, 14, dayLabel),
    S0: generateParams(0, 1, 0.01, 0.99, percentLabel),
    death_onset: generateParams(0, 1000, 1, 15, dayLabel),
    immunity_duration: generateParams(0.1, 10, 0.1, 1.0, yearLabel),
    life_expectancy: generateParams(1, 100, 1, 76, yearLabel),
    vaccination_rate: generateParams(0, 1, 0.01, 0.0, percentLabel)
};
