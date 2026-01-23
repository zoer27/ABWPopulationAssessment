#Code to fit the base model for ABW population assessment
#Base model includes JARPA, SOWER and photo-ID abundance estimates
#Low JARPA estimate is removed
#uses 2 * genetic Nmin for prior on N1973
#Sensitivities: 
    #code allows for multiple values of z
    #code allows for different priors on r
    #code allows for turning off data sources (JARPA, SOWER, photo-ID)
    #code allows for fixed survival rate to be included in model
#runs 4 chains, with random initial parameters
#can only do 1, 2, or 4 chains because each chain is hard coded in 
#does not do parallel sampling 
#uses slice sampling MCMC algorithm

#Zoe Rand
#last updated 8/15/25

using Optim
using Turing
using SliceSampling
using CSV
using DataFrames
using Distributions
using StatsPlots
using AdvancedMH
using MCMCChainsStorage
using HDF5
using Plots
using StatsPlots
#using TuringGLM #for negative binomial 2
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
years_Branch = abund_dat.Year[findall(==("Branch 2007"), abund_dat.Source)]
years_Ham = abund_dat.Year[findall(==("Hamabe 2023"), abund_dat.Source)]
#year indexes for abundance estimates
Br_index = findall(in(years_Branch), years)
Ham_index = findall(in(years_Ham), years)


#photo ID data
photo_est = CSV.read("data/photoID_abund.csv", DataFrame)
photo_idx = findall(in(2018), years)


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
FitThetaLogistic = function(r, r_star,lnN1973, catches, nyears, Br_index, Ham_index, z; fixed_s = false)
    #Input parameters (r, N1973, q_Ham)
    #N1973 = exp(lnN1973)
    #r = r
    #Find K using optimization
    FindK = function(K,r = r, r_star = r_star, N1973 = exp(lnN1973), catches = catches, z = z; fixed_s = fixed_s)
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

        if(fixed_s)
            surv = 0.96
             for i in 2:(1973-1904 +1)
                N[i] = N[i-1] + (1-surv)*r_star*N[i-1]*(1-(N[i-1]/K)^z) - catches[i-1]
                N[i] = posfun(N[i], 1)
                N[i] = max(N[i], 0) #shouldn't modify it, just for Julia to know it's positive
            end
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
    photo_n, photo_cv, photo_idx,z, niter, nchains, nthin, init_pars; 
    photo_id_switch = true, JARPA_switch = true, SOWER_switch = true, 
    prior_r = "uniform", fixed_s = false)
    #options for r prior are "uniform", "informative", or "calldensity"
   
    #Turing model
    @model function ThetaLogisticModel(N, CV, catches, nyears, Br_index, Ham_index,
       photo_n, photo_cv, photo_idx, z, 
       photo_id_switch, JARPA_switch, SOWER_switch, prior_r, fixed_s)
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
    
        #for posterior predictive:
        if(N === missing)
            N = Vector{Real}(undef, 12)
        end
        
        if(fixed_s)
            surv = 0.96 #fixed survival rate
            r_star = r/(1-surv)
        else
            r_star = r
        end
    
        
        #Get predictions
        K, Npred, NpredBranch, NpredHam = FitThetaLogistic(r, r_star, lnN1973, catches, nyears,
            Br_index, Ham_index, z; fixed_s = fixed_s)
        
       
        #predictions for photo id
        if(photo_id_switch)
            Nphoto = Npred[photo_idx]
            sigmaPhoto = sqrt.(log.((photo_cv.^2 .+ 1)))
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

        #generated quantities
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
        photo_n, photo_cv, photo_idx, z, 
        photo_id_switch, JARPA_switch, SOWER_switch, prior_r, fixed_s)
    
    #Sample from model
    if(nchains == 1)
        chains = sample(mod,externalsampler(RandPermGibbs(SliceSteppingOut(2.))),niter, 
            initial_params = init_pars, thinning = nthin)
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
        Ham_index,  Vector{Union{Missing, Float64}}(undef, length(photo_n)), photo_cv,
        photo_idx, z, photo_id_switch, JARPA_switch, SOWER_switch, prior_r, fixed_s)
    postpred = predict(post_mod, chains, include_all = false)
    
    #return samples and posterior predictive
    return(Chains = full_chains, PostPred = postpred)
    #return(full_chains)
end


pen = 0.0 #just in case it changed while testing

#All data
#setting up initial values
Random.seed!(123)
#right now it only works if num_values is 1, 2, or 4
num_values = 4  # Number of unique values you want to generate = number of chains
r_dist = Uniform(0.03, 0.09)  
r_inits = sample(rand(r_dist, 10*num_values), num_values; replace=false)

lnN1973_dist = Uniform(log(200), log(500))
ln1973_inits = sample(rand(lnN1973_dist, 10*num_values), num_values; replace=false)

#testing to see if they are reasonable
for i in 1:num_values
    K, N = FitThetaLogistic(r_inits[i],r_inits[i], ln1973_inits[i], catches, nyears, Br_index, Ham_index, 2.39)
    println(K)
end

add_cv_dist = Exponential(1)
add_cv_inits = sample(rand(add_cv_dist, 10*num_values), num_values; replace=false)

q_Ham_dist = Uniform(0, 0.5)
q_Ham_inits = sample(rand(q_Ham_dist, 10*num_values), num_values; replace=false)

init_pars = [[r_inits[i], ln1973_inits[i], add_cv_inits[i], add_cv_inits[i], q_Ham_inits[i]] for i in 1:num_values]

#number of iterations (after thinning)
niter =  1000
#thinning rate
nthinning = 50

#NOTE model names/numbers do not necessarily correlate to those in the manuscript

Telegram.sendMessage(tg, text = "Model 1 start")
 fit1 = MCMCFit(abund_dat.N, abund_dat.CV, catches, nyears, Br_index, Ham_index, 
     photo_est.N[2], photo_est.CV[2],  photo_idx, 2.39,  niter, num_values, nthinning, init_pars;
     photo_id_switch = true, JARPA_switch = true, SOWER_switch = true, prior_r = "uniform")

 h5open("Code/Results/Fit_BaseModel_uniformR_239.h5", "w") do f
     write(f, fit1.Chains)
 end

 #saving posterior predictive
 h5open("Code/Results/PostPred_BaseModel_uniformR_239.h5", "w") do f
     write(f, fit1.PostPred)
 end

 Telegram.sendMessage(tg, text = "Model 1 Complete")



Telegram.sendMessage(tg, text = "Model 4 Start")
 fit4 = MCMCFit(abund_dat.N, abund_dat.CV, catches, nyears, Br_index, Ham_index, 
     photo_est.N[2], photo_est.CV[2], photo_idx, 11.2,  niter, num_values, nthinning, init_pars;
     photo_id_switch = true, JARPA_switch = true, SOWER_switch = true, prior_r = "uniform")

h5open("Code/Results/Fit_BaseModel_uniformR_112.h5", "w") do f
     write(f, fit4.Chains)
 end

# #saving posterior predictive
 h5open("Code/Results/PostPred_BaseModel_uniformR_112.h5", "w") do f
     write(f, fit4.PostPred)
 end

 Telegram.sendMessage(tg, text = "Model 4 Complete")


Telegram.sendMessage(tg, text = "Model 6 Start")
 fit6 = MCMCFit(abund_dat.N, abund_dat.CV, catches, nyears, Br_index, Ham_index, 
     photo_est.N[2], photo_est.CV[2], photo_idx, 11.2,  niter, num_values, nthinning, init_pars;
     photo_id_switch = true, JARPA_switch = true, SOWER_switch = true, prior_r = "calldensity")

h5open("Code/Results/Fit_BaseModel_calldensityR_112.h5", "w") do f
     write(f, fit6.Chains)
 end

# #saving posterior predictive
 h5open("Code/Results/PostPred_BaseModel_calldensityR_112.h5", "w") do f
     write(f, fit6.PostPred)
 end

 Telegram.sendMessage(tg, text = "Model 6 Complete")



Telegram.sendMessage(tg, text = "Model 7 Start")
fit7 = MCMCFit(abund_dat.N, abund_dat.CV, catches, nyears, Br_index, Ham_index, 
    photo_est.N[2], photo_est.CV[2], photo_idx, 11.2,  niter, num_values, nthinning, init_pars;
    photo_id_switch = false, JARPA_switch = true, SOWER_switch = true, prior_r = "uniform")

h5open("Code/Results/Fit_BaseModel_nophoto.h5", "w") do f
    write(f, fit7.Chains)
end

#saving posterior predictive
h5open("Code/Results/PostPred_BaseModel_nophoto.h5", "w") do f
    write(f, fit7.PostPred)
end

Telegram.sendMessage(tg, text = "Model 7 Complete")



Telegram.sendMessage(tg, text = "Model 8 Start")
fit8 = MCMCFit(abund_dat.N, abund_dat.CV, catches, nyears, Br_index, Ham_index, 
    photo_est.N[2], photo_est.CV[2], photo_idx, 11.2,  niter, num_values, nthinning, init_pars;
    photo_id_switch = true, JARPA_switch = true, SOWER_switch = true, prior_r = "uniform")

h5open("Code/Results/Fit_BaseModel_noJARPA.h5", "w") do f
    write(f, fit8.Chains)
end

#saving posterior predictive
h5open("Code/Results/PostPred_BaseModel_noJARPA.h5", "w") do f
    write(f, fit8.PostPred)
end

Telegram.sendMessage(tg, text = "Model 8 Complete")


Telegram.sendMessage(tg, text = "Model 9 Start")
fit9 = MCMCFit(abund_dat.N, abund_dat.CV, catches, nyears, Br_index, Ham_index, 
   photo_est.N[2], photo_est.CV[2], photo_idx, 11.2,  niter, num_values, nthinning, init_pars;
    photo_id_switch = true, JARPA_switch = true, SOWER_switch = false, prior_r = "uniform")

h5open("Code/Results/Fit_BaseModel_noSOWER.h5", "w") do f
    write(f, fit9.Chains)
end 

#saving posterior predictive
h5open("Code/Results/PostPred_BaseModel_noSOWER.h5", "w") do f
    write(f, fit9.PostPred)
end

Telegram.sendMessage(tg, text = "Model 9 Complete")


Telegram.sendMessage(tg, text = "Model 10 Start")
fit10 = MCMCFit(abund_dat.N, abund_dat.CV, catches, nyears, Br_index, Ham_index, 
    photo_est.N[2], photo_est.CV[2], photo_idx, 11.2,  niter, num_values, nthinning, init_pars;
    photo_id_switch = false, JARPA_switch = false, SOWER_switch = true, prior_r = "uniform")


h5open("Code/Results/Fit_BaseModel_justSOWER.h5", "w") do f
    write(f, fit10.Chains)
end

#saving posterior predictive
h5open("Code/Results/PostPred_BaseModel_justSOWER.h5", "w") do f
    write(f, fit10.PostPred)
end

Telegram.sendMessage(tg, text = "Model 10 Complete")

#plot(fit10.Chains[:, [:r, :lnN1973, :q_Ham, :addCV_Branch, :addCV_Ham], :], )

Telegram.sendMessage(tg, text = "Model 11 Start")
fit11 = MCMCFit(abund_dat.N, abund_dat.CV, catches, nyears, Br_index, Ham_index, 
    photo_est.N[1], photo_est.CV[1], photo_idx, 11.2,  niter, num_values, nthinning, init_pars;
    photo_id_switch = true, JARPA_switch = true, SOWER_switch = true, prior_r = "uniform", fixed_s = true)

h5open("Code/Results/Fit_BaseModel_fixed_s.h5", "w") do f
    write(f, fit11.Chains)
end

#saving posterior predictive
h5open("Code/Results/PostPred_BaseModel_fixed_s.h5", "w") do f
    write(f, fit11.PostPred)
end

Telegram.sendMessage(tg, text = "Model 11 Complete")

plot(fit11.Chains[:, [:r, :lnN1973, :q_Ham, :addCV_Branch, :addCV_Ham], :], )

#For supplement: 
#Base model with matsuoka abundance estimates
abund_dat_mat = CSV.read("Data/Abund_est_mat.csv", DataFrame)
years_Branch = abund_dat_mat.Year[findall(==("Branch 2004"), abund_dat_mat.Source)]
years_Ham = abund_dat_mat.Year[findall(==("Matsuoka and Hakamada 2014"), abund_dat_mat.Source)]
#year indexes for abundance estimates
Br_index = findall(in(years_Branch), years)
Ham_index = findall(in(years_Ham), years)

Telegram.sendMessage(tg, text = "Model Mat Abund Start")
fit12 = MCMCFit(abund_dat_mat.N, abund_dat_mat.CV, catches, nyears, Br_index, Ham_index, 
    photo_est.N[2], photo_est.CV[2], photo_idx, 11.2,  niter, num_values, nthinning, init_pars;
    photo_id_switch = true, JARPA_switch = true, SOWER_switch = true, prior_r = "uniform")

h5open("Code/Results/Fit_BaseModel_uniformR_112_matabund.h5", "w") do f
    write(f, fit12.Chains)
end

#saving posterior predictive
h5open("Code/Results/PostPred_BaseModel_uniformR_112_matabund.h5", "w") do f
    write(f, fit12.PostPred)
end

Telegram.sendMessage(tg, text = "Model Mat Abund Complete")