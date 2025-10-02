// Parameter values for the SEIRS model
//
// This module generates discrete parameter value arrays for UI controls (sliders, dropdowns).
// Each parameter has a range of valid values with optional default and custom labels.

/**
 * Generates parameter configurations for model inputs
 * Creates an array of discrete values between start and end with specified step size.
 * Each value can be marked as default and formatted with a custom label.
 * 
 * @param {number} start - Starting value of the range
 * @param {number} end - Ending value of the range (inclusive)
 * @param {number} step - Increment step between consecutive values
 * @param {number|null} defaultVal - Value to mark as default (optional)
 * @param {Function|null} labelFunc - Optional function to format display labels
 * @returns {Array} Array of parameter objects with {value, default?, label?} structure
 */
const generateParams = (start, end, step, defaultVal = null, labelFunc = null) => {
    // Calculate number of steps needed to cover the range
    const numSteps = Math.round((end - start) / step);
    const params = [];
    
    // Generate each parameter value in the range
    for (let i = 0; i <= numSteps; i++) {
        const value = roundToTwoDecimals(start + i * step);
        const param = createParam(value, defaultVal, labelFunc);
        params.push(param);
    }
    
    return params;
};

/**
 * Rounds a number to two decimal places
 * Prevents floating-point precision errors in parameter generation
 * 
 * @param {number} num - Number to round
 * @returns {number} Number rounded to 2 decimal places
 */
const roundToTwoDecimals = (num) => Math.round(num * 100) / 100;

/**
 * Creates a parameter object with optional default flag and label
 * 
 * @param {number} value - The parameter value
 * @param {number|null} defaultVal - If matches value, marks this as default
 * @param {Function|null} labelFunc - Optional function to generate display label
 * @returns {Object} Parameter object with value and optional default/label properties
 */
const createParam = (value, defaultVal, labelFunc) => {
    const param = { value };
    
    // Mark as default if this value matches the default
    if (value === defaultVal) {
        param.default = true;
    }
    
    // Apply custom label formatter if provided
    if (labelFunc) {
        param.label = labelFunc(value);
    }
    
    return param;
};

// ============================================================================
// Label Formatters
// ============================================================================
// These functions format numerical values into human-readable labels for the UI

/**
 * Formats a value as days with proper pluralization
 * @example dayLabel(1) => "1 day", dayLabel(5) => "5 days"
 */
const dayLabel = (val) => `${val} day${val > 1 ? 's' : ''}`;

/**
 * Formats a value as years with proper pluralization
 * @example yearLabel(1) => "1 year", yearLabel(10) => "10 years"
 */
const yearLabel = (val) => `${val} year${val > 1 ? 's' : ''}`;

/**
 * Formats a proportion (0-1) as a percentage
 * @example percentLabel(0.5) => "50%", percentLabel(0.99) => "99%"
 */
const percentLabel = (val) => `${Math.round(val * 100)}%`;
const percentLabel2 = (val) => `${val}%`;

const logScaleLabel = (val) => `${val > 0 ? 'Log10' : 'Linear'}`;

// ============================================================================
// Parameter Value Definitions
// ============================================================================
// Defines the valid ranges and defaults for all SEIR model parameters

export const param_vals = {
    // Simulation duration: 50 to 3000 days in 50-day increments (default: 3000)
    n_days: generateParams(50, 3000, 50, 3000, dayLabel),
    
    // Y-axis maximum for plot: 5% to 100% in 5% increments (default: 100%)
    y_max: generateParams(5, 100, 5, 100, percentLabel2),
    
    // Basic reproduction number: 1.0 to 5.0 in 0.1 increments (default: 3.0)
    // R0 represents average number of secondary infections from one infected individual
    R0: generateParams(1, 5, 0.1, 3.0),
    
    // Latent period: 1 to 30 days (default: 7 days)
    // Time from exposure to becoming infectious
    latent_period: generateParams(1, 30, 1, 7, dayLabel),
    
    // Infectious period: 1 to 30 days (default: 14 days)
    // Duration of infectiousness before recovery
    infectious_period: generateParams(1, 30, 1, 14, dayLabel),
    
    // Initial susceptible proportion: 0% to 100% in 1% increments (default: 99%)
    S0: generateParams(0, 1, 0.01, 0.99, percentLabel),
    
    // Death onset: 0 to 1000 days (default: 0, meaning no disease-induced death)
    // Average time from infection to death for fatal cases
    death_onset: generateParams(0, 1000, 1, 0, dayLabel),
    
    // Immunity duration: 0.1 to 10 years in 0.1-year increments (default: 1 year)
    // How long recovered individuals remain immune before becoming susceptible again
    immunity_duration: generateParams(0, 10, 0.1, 1.0, yearLabel),
    
    // Life expectancy: 1 to 100 years (default: 76 years)
    // Average lifespan, used to calculate natural birth/death rate
    life_expectancy: generateParams(1, 100, 1, 76, yearLabel),
    
    // Vaccination rate: 0% to 100% in 1% increments (default: 0%)
    // Proportion of population vaccinated (enter recovered compartment)
    vaccination_rate: generateParams(0, 1, 0.01, 0.0, percentLabel),

    use_log_scale: generateParams(0, 1, 1, 0, logScaleLabel)
};