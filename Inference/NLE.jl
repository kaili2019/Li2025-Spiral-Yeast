# This script implement NLE for inference of the ABM of sprial yeast colony (M. Magnusii).
# 
#
# Kai Li
# 29 Jan 2025

using CSV, Tables
using DataFrames
using Plots
using StatsPlots

using StatsBase
using Flux, Distributions, MixtureDensityNetworks
using JLD2

logitAWRI50um = CSV.read("ss_galaxy_r3_23June2025.csv", DataFrame; delim = ',')

logit(x) = @. log(x/(1-x));
invlogit(x) = @. 1/(1+exp(-x));

# Programmatically create the formula
predictor_names = names(logitAWRI50um)[7:end]  

new_data = logitAWRI50um[:,predictor_names]

theta_names = ["nStar","p2sProb","s2pProb","gamma","pa","angle_prolif"];

Y = Matrix(new_data)

# apply inverse logit
X = invlogit(reshape([logitAWRI50um.nStar; logitAWRI50um.p2sProb; logitAWRI50um.s2pProb; logitAWRI50um.gamma; logitAWRI50um.pa; logitAWRI50um.angle_prolif],length(logitAWRI50um.nStar),6))

# train model
const epochs = 600
const batchsize = 16
const mixtures = 6
const layers = [64, 64]

model = MixtureDensityNetwork(6, 6, layers, mixtures);

model, report = MixtureDensityNetworks.fit!(model, Matrix((X')), Matrix((Y')); epochs=epochs, opt=Flux.Adam(1e-3), batchsize=batchsize)

galaxy_model = Dict()
galaxy_model["model"] = model
galaxy_model["report"] = report
# Save the entire dict to a file
jldsave("galaxy_models_23June2025_invlogit_r3.jld2"; galaxy_model)

plot(1:epochs, report.learning_curve, xlabel="Epochs", ylabel="Loss")

sel_idx = 1;

ground_truth_50um = CSV.read("ss_galaxy_exp_r3_ss_23June2025.csv", DataFrame; delim = ',')

ground_truth = ground_truth_50um[sel_idx,predictor_names]

ss_mean_500um = (collect(ground_truth))

target_ss = ss_mean_500um;

function log_likelihood(sample)

    # take slice of likehood from MDN 
    cond = model(reshape(sample, (6,1)))[1] 
    # evaluate the summary statistics values at pdf 
    likelihood_val = logpdf(cond,target_ss)

    return likelihood_val

end

# Beta prior densities 

A = [1, 1, 1, 1, 1, 1]; 
B = [1, 1, 1, 1, 1, 1];

# sample from prior 
function log_prior(sample)

    Beta_prior = logpdf.(Beta.(A,B),reshape((sample),6,1));

    return sum(Beta_prior);

end

function log_post_sample(sample)

    if any(x -> x < 0 || x > 1, sample)
        return -Inf
    end

    log_likelihood(sample)+log_prior(sample);

end

# number of samples to take 
N = 200000;
# matrix to store samples of theta 
samples = Matrix{Float64}(undef,N,6); 
samples[1,:] = ([0.2, 0.5, 0.5, 0.5, 0.5, 0.5]);
store_alpha = Matrix{Float64}(undef,N,1);

store_proposal = Matrix{Float64}(undef,N,6); 
store_proposal[1,:] = samples[1,:];

# MCMC sampler 
for ii = 2:N 

    # sample from proposal 
    proposal_thetaStar = rand.(Normal.(samples[ii-1,:], 0.1));
    store_proposal[ii,:] = proposal_thetaStar;

    # calculate acceptance probability  
    numerator = log_post_sample(proposal_thetaStar);
    
    denominator = log_post_sample(samples[ii-1,:]);

    alpha = numerator[1]-denominator[1];

    store_alpha[ii-1] = alpha;
    
    # accept with probability alpha 
    if log(rand()) < alpha 
        samples[ii,:] = proposal_thetaStar;
    else
        samples[ii,:] = samples[ii-1,:];
    end

end

stephist((samples[10000:end,1:6]),alpha=1,normalize=:pdf)
mean(eachrow((samples)))

plot((samples[10000:200000,:]))

# overlay prior with accepted samples 
names_vec = ["nStar", "Pps", "Psp", "gamma", "pa", "angle_prolif"];

plots = [];
for idx = 1:6
    p = histogram((samples[10000:5:end,idx]),normalize=:pdf,label=names_vec[idx],linealpha=0,legend=:topright,bins=0:0.02:1)
    plot!(p, x -> pdf.(Beta.(A[idx],B[idx]), x),0:0.01:1,label="Prior",lw=3,legend=:topright) 
    scatter!(p, [mean(eachrow((samples)))[idx]], [0.2], label="MDN sample mean", color="red", marker=:star, ms=8,legend=:topright)

    if idx == 6
        p = histogram((samples[10000:5:end,idx])*pi/20,normalize=:pdf,label=names_vec[idx],linealpha=0,legend=:topright,bins=0:0.002:pi/20)
        plot!(p, x -> pdf.(Uniform.(0,pi/20), x), 0:0.01:pi/20,label="Prior",lw=3,legend=:topright) 
        scatter!(p, [mean(eachrow((samples)))[idx]]*pi/20, [0.2], label="MDN sample mean", color="red", marker=:star, ms=8,legend=:topright)
    end
    push!(plots, p)
end

# p1 = plot(plots...,layout=(2,3),size=(1000,600),suptitle="Prior Predictive Checking")

p1 = plot(plots...,layout=(2,3),size=(1000,600))

# savefig(p1, "PPC_Sample_23June2025_galaxy_r3_v2.pdf")

# CSV.write("samples_23Jun2025_galaxy_r3_v2.csv",Tables.table((samples[10000:5:end,:])),writeheader=false)
