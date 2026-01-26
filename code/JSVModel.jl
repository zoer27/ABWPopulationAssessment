#Code to fit the JSV model for ABW population assessment
#Same as base model but includes JSV data and likelihood
#Zoe Rand
#last updated 12/5/25
#Note that file paths for storing results should be specified

using Optim
using Turing
using AdvancedMH
using SliceSampling
using CSV
using DataFrames
using Distributions
using StatsPlots
using MCMCChainsStorage
using HDF5
using Plots
using StatsPlots
using TuringGLM #for negative binomial 2
using Random 

#setting up Telegram for notifications
using Telegram
bot_token = "7451303674:AAF3k_1zVp9wWdWnvO7q72nZncx7T_yP8K4"

chat_id = "7281464972"

tg = TelegramClient(bot_token, chat_id = chat_id)

#Read in Data
## abundance estimates
abund_dat = CSV.read("data/Abund_Est.csv", DataFrame)

#taking out low JARPA estimate 
abund_dat = abund_dat[abund_dat.N .> 6, :]

## catch data
catchdf = CSV.read("data/ABW_catch_series_2024.csv", DataFrame)
catches = deleteat!(catchdf.ABW_catch, 1) #removing first row because it's 1903

#years
years = (minimum(catchdf.Year)+1):2024
nyears = length(years)

#adding 0s to catches
catches = [catches; zeros(nyears - length(catches))]

#year indexes for abundance estimates
    # Get the Year entries from abund_dat when Source = "Branch 2004"
years_Branch = abund_dat.Year[findall(==("Branch 2004"), abund_dat.Source)]
years_Ham = abund_dat.Year[findall(==("Hamabe 2023"), abund_dat.Source)]
#year indexes for abundance estimates
Br_index = findall(in(years_Branch), years)
Ham_index = findall(in(years_Ham), years)


#photo ID data
photo_est = CSV.read("data/PhotoID_abund.csv", DataFrame)
photo_idx = findall(in(2018), years)


#which survival rate to use
photodf_run = photodf[findall(in(0.96),photodf.Survival ), :] 


#JSV data
JSVdf = CSV.read("data/JSV.csv", DataFrame)
sightings_early = JSVdf.Sightings[JSVdf.Year .< 1978]
sightings_late = JSVdf.Sightings[JSVdf.Year .>= 1978]
effort_early = JSVdf.Effort[JSVdf.Year .< 1978]
effort_late = JSVdf.Effort[JSVdf.Year .>= 1978]
#dividing effort by 1000 to see if this helps with numerical stability
effort_early = effort_early ./ 100000
effort_late = effort_late ./ 100000

JSV_early_idx = findall(in(JSVdf.Year[JSVdf.Year .< 1978]), years)
JSV_late_idx = findall(in(JSVdf.Year[JSVdf.Year .>= 1978]), years)


#posfun function
posfun = function(x, eps)
    global pen
    if ( x >= eps )
        return x
    else 
        pen += 0.1 * (x-eps)^2
        return (eps/(2.0-x/eps))
    end
end

#global variable for penalty
pen = 0.0
#Function to fit the model
FitThetaLogistic = function(r, lnN1973, catches, nyears, Br_index, Ham_index, z)
    #Input parameters (r, N1973, q_Ham)
    #N1973 = exp(lnN1973)
    #r = r
    #Find K using optimization
    FindK = function(K,r = r, N1973 = exp(lnN1973), catches = catches, z = z)
        ##Given catches, and r,  and assuming 1904 = proposed K,
        #### calculate population trajectory from 1904 to 1973 and compare to
        ##### 1973--calculate the difference squared
        global pen
        pen = 0.0 #start penalty at 0
        N = Vector{Any}(undef, 1973-1904 +1)
        N[1] = K
        for i in 2:(1973-1904 +1)
            N[i] = N[i-1] + r*N[i-1]*(1-(N[i-1]/K)^z) - catches[i-1]
            N[i] = posfun(N[i], 1)
            N[i] = max(N[i], 0) #shouldn't modify it, just for Julia to know it's positive
        end
        #println(N)
       #println(pen)
    return((N[70]-N1973)^2 + pen)
    end
    res = optimize(FindK, 100000, 400000, abs_tol = 10^-20)
    #print(Optim.abs_tol(res))
    K = res.minimizer
    #print(K)
    ##Using r and K from minimization, calculate full population trajectory
    N = Vector{Any}(undef, nyears)
    N[1] = K
    for i in 2:nyears
        N[i] = N[i-1] + r*N[i-1]*(1-(N[i-1]/K)^z) - catches[i-1]
        N[i] = max(N[i], 0.001)
    end
    ## Return predictions  that correspond to data
    return(K, N, N[Br_index], N[Ham_index])
    #return(K, N, N[Br_index], N[Ham_index], N[70]) #add in N[70] which is N1973
end



#test above
#FitThetaLogistic(0.07, log(100.0), catches, nyears, Br_index, Ham_index, 2.39)

#MCMC Function
MCMCFit = function(N, CV, catches, nyears, Br_index, Ham_index,
    sightings_early, sightings_late,
    effort_early, effort_late, JSV_early_idx, JSV_late_idx,
    photo_n, photo_cv, photo_idx,z, nsamp, nchains, nthin, init_pars; 
    photo_id_switch = true, JARPA_switch = true, SOWER_switch = true, JSV_switch = true,
    prior_r = "uniform")
    #options for r prior are "uniform", "informative", or "calldensity"
   
    #Turing model
    @model function ThetaLogisticModel(N, CV, catches, nyears, Br_index, Ham_index,
        sightings_early, sightings_late,
        effort_early, effort_late, JSV_early_idx, JSV_late_idx,
        photo_n, photo_cv, photo_idx, z, 
        photo_id_switch, JARPA_switch, SOWER_switch, JSV_switch, prior_r)
        #Priors
        if(prior_r == "uniform")
            r ~ Uniform(0, 0.114)
        elseif(prior_r == "informative")
            r ~ Normal(0.043, 0.019) #from Branch 2004
        elseif(prior_r == "calldensity")
            #from call density analysis
            #annual mean increase was 3.73
            #monthly standard deviation was 0.00428
            #annual sd is sqrt(12)*0.00428 = 0.0148
            r ~ Normal(0.0373, 0.0148)
        end

        #prior based on 2*nMin from mtDNA
        lnN1973 ~ Uniform(log(106), log(2000))

        if(JARPA_switch)
            q_Ham ~ Uniform(0, 1)
        end

        if(SOWER_switch)
            addCV_Branch ~ Exponential(1)
        end
        if(JARPA_switch)
            addCV_Ham ~ Exponential(1)
        end
        
        if(JSV_switch)
            theta_early ~ Gamma(0.01, 10)
            theta_late ~ Gamma(0.01, 10)
            q_early ~ Uniform(0,1)
            q_late ~ Uniform(0,1)
        end

        #for posterior predictive:
        if(N === missing)
            N = Vector{Real}(undef, 12)
        end
        
        #Get predictions
        K, Npred, NpredBranch, NpredHam = FitThetaLogistic(r, lnN1973, catches, nyears, 
            Br_index, Ham_index, z)
        
       
        #predictions for photo id
        if(photo_id_switch)
            Nphoto = Npred[photo_idx]
            sigmaPhoto = sqrt.(log.((photo_cv.^2 .+ 1)))
        end

        #predictions for JSV
        if(JSV_switch)
            NJSV_early = Npred[JSV_early_idx]
            NJSV_late = Npred[JSV_late_idx]
            Lambda_early = (NJSV_early.*effort_early .* q_early)
            Lambda_late = (NJSV_late.*effort_late .* q_late)
        end

        #Calculate sigma for abundance
        CV_squared = CV.^2
        CV_squared_Br = CV_squared[1:3]
    
        CV_squared_Ham = CV_squared[4:end]
        
        if(SOWER_switch)
            ADDCV_Br_squared = addCV_Branch^2
            sigmaBranch = sqrt.(log.(CV_squared_Br .+ ADDCV_Br_squared .+ 1))
        end

        if(JARPA_switch)
            AddCV_Ham_squared = addCV_Ham^2
            sigmaHam = sqrt.(log.(CV_squared_Ham .+ AddCV_Ham_squared .+ 1))
        end
            
        
        #Likelihood
        #SOWER
        if(SOWER_switch)
            for i in eachindex(Br_index)
                N[i] ~ LogNormal(log(NpredBranch[i]), sigmaBranch[i])
            end
        end
        #JARPA

        if(JARPA_switch)
            for i in 3 .+ eachindex(Ham_index)
                N[i] ~ LogNormal(log(q_Ham*NpredHam[i-3]), sigmaHam[i-3])
            end
        end
    

        #photoID
        if(photo_id_switch)
            for i in eachindex(photo_n)
                photo_n[i] ~ LogNormal(log(Nphoto[i]), sigmaPhoto[i])
            end
        end

        #JSV
        if(JSV_switch)
            for i in eachindex(sightings_early)
               sightings_early[i] ~ TuringGLM.NegativeBinomial2(Lambda_early[i], theta_early)
            end
            for i in eachindex(sightings_late)
                sightings_late[i] ~ TuringGLM.NegativeBinomial2(Lambda_late[i], theta_late)
            end
        end

        Npopkeys = ["Npop$(i)" for i in 1:nyears] #creating named tuple for Npop for generated quantities
        Npred_namedtuple = NamedTuple{Tuple(Symbol.(Npopkeys))}(Tuple(Npred))
        Nmin = minimum(Npred)
        Nmin_div_K = minimum(Npred)/K
        N2024 = Npred[end]
        N2024_div_K = N2024/K
        return(K = K, Nmin = Nmin, Nmin_div_K = Nmin_div_K, 
        N2024 = N2024, N2024_div_K = N2024_div_K, Npred_namedtuple...) #generated quantites
    end

    #Declare model
    mod = ThetaLogisticModel(N, CV, catches, nyears, Br_index, Ham_index,
            sightings_early, sightings_late,
            effort_early, effort_late, JSV_early_idx, JSV_late_idx,
            photo_n, photo_cv, photo_idx, z, 
            photo_id_switch, JARPA_switch, SOWER_switch, JSV_switch, prior_r)
    
    #Sample from model
    if(nchains == 1)
        chains = sample(mod,externalsampler(RandPermGibbs(SliceSteppingOut(2.))), niter, 
            initial_params = init_pars[1], thinning = nthin)
    end
    if(nchains == 2)
        chain1 = sample(mod,externalsampler(RandPermGibbs(SliceSteppingOut(2.))),niter, 
            initial_params = init_pars[1], thinning = nthin)
        chain2 = sample(mod,externalsampler(RandPermGibbs(SliceSteppingOut(2.))),niter,
            initial_params = init_pars[2], thinning = nthin)
        chains = chainscat(chain1, chain2)
    end

    if(nchains == 4)
        chain1 = sample(mod,externalsampler(RandPermGibbs(SliceSteppingOut(2.))),niter, 
            initial_params = init_pars[1], thinning = nthin)
        chain2 = sample(mod,externalsampler(RandPermGibbs(SliceSteppingOut(2.))),niter,
            initial_params = init_pars[2], thinning = nthin)
        chain3 = sample(mod,externalsampler(RandPermGibbs(SliceSteppingOut(2.))),niter, 
            initial_params = init_pars[3], thinning = nthin)
        chain4 = sample(mod,externalsampler(RandPermGibbs(SliceSteppingOut(2.))),niter,
            initial_params = init_pars[4], thinning = nthin)

        chains = chainscat(chain1, chain2, chain3, chain4)
    end
    #Add generated quantities
    full_chains = append_generated_quantities(chains, Turing.generated_quantities(mod, chains))
    
    #get posterior predictive
    post_mod = ThetaLogisticModel(missing, CV, catches, nyears, Br_index, 
        Ham_index, Vector{Union{Missing, Float64}}(undef, length(sightings_early)), Vector{Union{Missing, Float64}}(undef, length(sightings_late)),
        effort_early, effort_late, JSV_early_idx, JSV_late_idx,  Vector{Union{Missing, Float64}}(undef, length(photo_n)), photo_cv,
        photo_idx, z, photo_id_switch, JARPA_switch, SOWER_switch, JSV_switch, prior_r)
    postpred = predict(post_mod, chains, include_all = false)
    
    #return samples and posterior predictive
    return(Chains = full_chains, PostPred = postpred)
    #return(full_chains)
end

#Threads.nthreads()
pen = 0.0 #just in case it changed while testing

#All data
#setting up initial values
Random.seed!(445)
num_values = 4  # Number of unique values you want to generate = number of chains, currently only works for 1, 2, or 4
r_dist = Uniform(0.03, 0.09)  
r_inits = sample(rand(r_dist, 10*num_values), num_values; replace=false)

lnN1973_dist = Uniform(log(200), log(500))
ln1973_inits = sample(rand(lnN1973_dist, 10*num_values), num_values; replace=false)

#testing to see if they are reasonable
for i in 1:num_values
    K, N = FitThetaLogistic(r_inits[i], ln1973_inits[i], catches, nyears, Br_index, Ham_index, 2.39)
    println(K)
end

add_cv_dist = Exponential(1)
add_cv_inits = sample(rand(add_cv_dist, 10*num_values), num_values; replace=false)

q_Ham_dist = Uniform(0, 0.5)
q_Ham_inits = sample(rand(q_Ham_dist, 10*num_values), num_values; replace=false)

theta_dist = Gamma(0.01, 10)
theta_early_inits = sample(rand(theta_dist, 10*num_values), num_values; replace=false) 
theta_late_inits = sample(rand(theta_dist, 10*num_values), num_values; replace=false)
q_JSV_dist = Uniform(0, 1)
q_JSV_early_inits = sample(rand(q_JSV_dist, 10*num_values), num_values; replace=false) 
q_JSV_late_inits = sample(rand(q_JSV_dist, 10*num_values), num_values; replace=false)

init_pars = [[r_inits[i], ln1973_inits[i], add_cv_inits[i], add_cv_inits[i], q_Ham_inits[i],
            theta_early_inits[i], theta_late_inits[i], q_JSV_early_inits[i], q_JSV_late_inits[i]] for i in 1:num_values]

#number of iterations (after thinning)
niter = 1000
#thinning rate
nthinning = 50


Telegram.sendMessage(tg, text = "Model JSV Start")
fit4 = MCMCFit(abund_dat.N, abund_dat.CV, catches, nyears, Br_index, Ham_index, 
    sightings_early, sightings_late, effort_early, effort_late, JSV_early_idx, JSV_late_idx,
   photo_est.N[2], photo_est.CV[2], photo_idx, 11.2,  niter, num_values, nthinning, init_pars;
    photo_id_switch = true, JARPA_switch = true, SOWER_switch = true, JSV_switch = true, prior_r = "uniform")

h5open("Code/Results/Fit_JSVModel_uniformR_112.h5", "w") do f
    write(f, fit4.Chains)
end

#saving posterior predictive
h5open("Code/Results/PostPred_JSVModel_uniformR_112.h5", "w") do f
    write(f, fit4.PostPred)
end

Telegram.sendMessage(tg, text = "Model JSV Complete")

