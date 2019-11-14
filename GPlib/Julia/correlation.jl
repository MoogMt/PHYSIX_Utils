module correlation

using Statistics

function autocorrAvg( signal::Vector{T1}, max_lag::T2 ) where { T1 <: Real, T2 <: Int }
    nb_step=size(signal)[1]
    autocor_signal=zeros(max_lag)
    avg_sig=mean(signal)
    # Loop over tau
    for lag=0:max_lag-1
        for step=1:nb_step-lag
            autocor_signal[lag+1] += (signal[step]-avg_sig)*(signal[step+lag]-avg_sig)
        end
        autocor_signal[lag+1] /= (max_lag-lag)
    end
    return autocor_signal
end

function autocorrAvgSig( signal::Vector{T1}, max_lag::T2 ) where { T1 <: Real, T2 <: Int }
    nb_step=size(signal)[1]
    autocor_signal=zeros(max_lag)
    avg_sig=mean(signal)
    var_sig=Statistics.var(signal)
    # Loop over tau
    for lag=0:max_lag-1
        for step=1:nb_step-lag
            autocor_signal[lag+1] += (signal[step]-avg_sig)*(signal[step+lag]-avg_sig)
        end
        autocor_signal[lag+1] /= ((max_lag-lag)*var_sig)
    end
    return autocor_signal
end

function autocorr( signal::Vector{T1}, max_lag::T2 ) where { T1 <: Real, T2 <: Int }
    nb_step=size(signal)[1]
    autocor_signal=zeros(max_lag)
    # Loop over tau
    for lag=0:max_lag-1
        for step=1:nb_step-lag
            autocor_signal[lag+1] += (signal[step])*(signal[step+lag])
        end
        autocor_signal[lag+1] /= (max_lag-lag)
    end
    return autocor_signal
end

function autocorrNorm( signal::Vector{T1}, max_lag::T2 ) where { T1 <: Real, T2 <: Int }
    autocor_signal = autocorr(signal,max_lag)
    autocor_signal/=autocor_signal[1]
    return autocor_signal
end

end
