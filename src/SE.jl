using FFTW
using Polynomials
using LinearAlgebra
using Statistics

function detrend(x,y; deg = 1)
    fit = polyfit(x,y,deg)
    detrended = y - fit(x)
    return detrended
end

function smooth(data; m = 5)
    if m != 0
        smoothed_data = Float64[]
        weigths = Int[]
        for w in collect(-m:m)
            p = (abs(w)+1)*2
            if p <= 20
                append!(weigths,p)
            else
                append!(weigths,20)
            end
        end
        for i in 1:length(data)
            if i > m && i < length(data) - m
                smoothed_Iw = 0
                for (l,k) in enumerate(collect(-m:m))
                    smoothed_Iw += (1/weigths[l])*data[i+k]
                end
                append!(smoothed_data, smoothed_Iw)
            end
        end
        return smoothed_data
    elseif m == 0
        return data
    end
end

"""
    vectorize(data)
    
Realizes the one-hot encoding of the provided time-serie.
Each category is assigned to it's own vector. for example,
if we provided [1,2,3,2,1,2,3] as input data, we get 
[[1,0,0],[0,1,0],[0,0,0],[1,0,0],,[0,1,0],[0,0,0]] as output.
"""

function vectorize(data)
    categories = unique(data)
    deleteat!(categories,length(categories))
    sorted_data = zeros(length(categories),length(data))
    for (t,d) in enumerate(data)
        for (i,j) in enumerate(categories)
            if d == j
                sorted_data[i,t] = 1
            end
        end
    end
    return sorted_data
end

function periodogram_matrix(ts;smoothing_degree=3)
    periodo = Array{Float64,3}(undef,length(ts[:,1]),length(ts[:,1]),length(ts[1,:]))
    for i in 1:length(ts[:,1])
        for j in 1:length(ts[:,1])
            tmp = smooth(real((1/length(ts[1,:]))*fft(ts[i,:]).*(conj(fft(ts[j,:])))); m= smoothing_degree)
            for w in 1:length(ts[1,:])-(2*smoothing_degree+1)
                if tmp[w] == NaN || tmp[w] == Inf
                    print("Unexpected NaN or Inf at position", w)
                else
                periodo[i,j,w] = tmp[w]
                end
            end
        end
    end

    return periodo
end

function varcov(ts)
    cov_matrix = zeros(length(ts[:,1]),length(ts[:,1]))
    for i in 1:length(ts[:,1])
        for j in 1:length(ts[:,1])
            if cov_matrix[i,j] == NaN || cov_matrix[i,j] == Inf
                print("Unexpected NaN or Inf at position", i,j)
            else
                cov_matrix[i,j] = sum((1/(length(ts[1,:]).-1))*(ts[i,:].-mean(ts[i,:])).*(ts[j,:].-mean(ts[j,:])))
            end
        end
    end
    return cov_matrix
end


function spectral_function(ts; m = 3)
    result = []
    eigvec_any = []
        S = sqrt(varcov(vectorize(ts)))
        p = periodogram_matrix(vectorize(ts); smoothing_degree = m)
        for i in 1:trunc(Int,length(ts)/2)
            if any(isnan.(S*p[:,:,i]*S))
                print("unexpected NaN at position : ", i,"\n", S,"\t")
            else
                append!(result,findmax(real(eigvals(S*p[:,:,i]*S)))[1])
                append!(eigvec_any,real(eigvecs(S*p[:,:,i]*S))[findmax(real(eigvals(S*p[:,:,i]*S)))[2],:])
            end
        end
        eigvec = reshape(Array{Float64}(eigvec_any),length(unique(ts))-1,:)'
        return collect(1:length(result))/(2*length(result)),result,eigvec
end

function scaling(ts;start = 1, stop = 500)
    dat = spectral_function(ts)
    pos = findmax(ep[2][start:stop])[2]
    return dat[3][pos,:]
end

function findmax_in(xserie,yserie,xlim)
    range = []
    real_pos = []
    for (index,value) in enumerate(xserie)
        if value>= xlim[1] && value <= xlim[2]
            append!(range,yserie[index])
            append!(real_pos,index)
        end
    end
    max,pos = findmax(range)
    return max, xserie[real_pos[pos]], real_pos[pos]
end

export spectral_function, findmax_in, detrend, smooth
