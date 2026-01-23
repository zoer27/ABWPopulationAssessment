#This code reads in the saved model fits generated from the 
#"BaseModelandSensitivities.jl" as well as "JSVModel.jl"
#code and checks diagnostics and saves the posteriors for plotting
#and results

using Turing
using CSV
using DataFrames
using Distributions
using StatsPlots
using AdvancedMH
using MCMCChainsStorage
using HDF5
using Plots
using StatsPlots
using MCMCChains

#read in fit
fit1 = h5open("Code/Results/Fit_BaseModel_uniformR_239.h5", "r") do f
    read(f, Chains)
end
out1 = DataFrame(fit1[501:end, :, :])

#remove burnin and save
CSV.write("Code/Results/Fit_BaseModel_uniformR_239.csv", out1)


#read in posterior predictiver
postpred1 = h5open("Code/Results/PostPred_BaseModel_uniformR_239.h5", "r") do f
    read(f, Chains)
end
postpredout1 = DataFrame(postpred1[501:end, :, :])

#save posterior predictive
CSV.write("Code/Results/PostPred_BaseModel_uniformR_239.csv", postpredout1)


#read in fit
fit3 = h5open("Code/Results/Fit_BaseModel_calldensityR_239.h5", "r") do f
    read(f, Chains)
end

plot(fit3[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham"],:])
autocorplot(fit3[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham"],:])
out3 = DataFrame(fit3[501:end, :, :])

#remove burnin and save
CSV.write("Code/Results/Fit_BaseModel_calldensityR_239.csv", out3)

#read in posterior predictiver
postpred3 = h5open("Code/Results/PostPred_BaseModel_calldensityR_239.h5", "r") do f
    read(f, Chains)
end
postpredout3 = DataFrame(postpred3[501:end, :, :])

#save posterior predictive
CSV.write("Code/Results/PostPred_BaseModel_calldensityR_239.csv", postpredout3)

#read in fit
fit4 = h5open("Code/Results/Fit_BaseModel_uniformR_112.h5", "r") do f
    read(f, Chains)
end

plot(fit4[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham"],:])
autocorplot(fit4[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham"],:])
out4 = DataFrame(fit4[501:end, :, :])

#remove burnin and save
CSV.write("Code/Results/Fit_BaseModel_uniformR_112.csv", out4)

#read in posterior predictiver
postpred4 = h5open("Code/Results/PostPred_BaseModel_uniformR_112.h5", "r") do f
    read(f, Chains)
end

postpredout4 = DataFrame(postpred4[501:end, :, :])

#save posterior predictive
CSV.write("Code/Results/PostPred_BaseModel_uniformR_112.csv", postpredout4)


#read in fit
fit6 = h5open("Code/Results/Fit_BaseModel_calldensityR_112.h5", "r") do f
    read(f, Chains)
end

plot(fit6[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham"],:])
autocorplot(fit6[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham"],:])
out6 = DataFrame(fit6[501:end, :, :])

#remove burnin and save
CSV.write("Code/Results/Fit_BaseModel_calldensityR_112.csv", out6)

#read in posterior predictiver
postpred6 = h5open("Code/Results/PostPred_BaseModel_calldensityR_112.h5", "r") do f
    read(f, Chains)
end
postpredout6 = DataFrame(postpred6[501:end, :, :])

#save posterior predictive
CSV.write("Code/Results/PostPred_BaseModel_calldensityR_112.csv", postpredout6)


#read in fit
fit10 = h5open("Code/Results/Fit_JSVModel_uniformR_112.h5", "r") do f
    read(f, Chains)
end

plot(fit10[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham", "theta_early", "theta_late", "q_early", "q_late"],:])
autocorplot(fit10[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham", "theta_early", "theta_late", "q_early", "q_late"],:])
out10 = DataFrame(fit10[501:end, :, :])

#remove burnin and save
CSV.write("Code/Results/Fit_JSVModel_uniformR_112.csv", out10)
#read in posterior predictiver
postpred10 = h5open("Code/Results/PostPred_JSVModel_uniformR_112.h5", "r") do f
    read(f, Chains)
end
postpredout10 = DataFrame(postpred10[501:end, :, :])
#save posterior predictive
CSV.write("Code/Results/PostPred_JSVModel_uniformR_112.csv", postpredout10)



fit13 = h5open("Code/Results/Fit_BaseModel_noSOWER.h5", "r") do f
    read(f, Chains)
end
plot(fit13[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Ham"],:])
autocorplot(fit13[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Ham"],:])
out13 = DataFrame(fit13[501:end, :, :])
#remove burnin and save
CSV.write("Code/Results/Fit_BaseModel_noSOWER.csv", out13)

#read in posterior predictiver
postpred13 = h5open("Code/Results/PostPred_BaseModel_noSOWER.h5", "r") do f
    read(f, Chains)
end
postpredout13 = DataFrame(postpred13[501:end, :, :])    
#save posterior predictive
CSV.write("Code/Results/PostPred_BaseModel_noSOWER.csv", postpredout13)

fit14 = h5open("Code/Results/Fit_BaseModel_noJARPA.h5", "r") do f
    read(f, Chains)
end

plot(fit14[501:end, ["r", "lnN1973", "K", "addCV_Branch"],:])
autocorplot(fit14[501:end, ["r", "lnN1973", "K", "addCV_Branch"],:])
out14 = DataFrame(fit14[501:end, :, :])
#remove burnin and save
CSV.write("Code/Results/Fit_BaseModel_noJARPA.csv", out14)  

#read in posterior predictiver
postpred14 = h5open("Code/Results/PostPred_BaseModel_noJARPA.h5", "r") do f
    read(f, Chains)
end
postpredout14 = DataFrame(postpred14[501:end, :, :])
#save posterior predictive
CSV.write("Code/Results/PostPred_BaseModel_noJARPA.csv", postpredout14)

fit15 = h5open("Code/Results/Fit_BaseModel_nophoto.h5", "r") do f
    read(f, Chains)
end
plot(fit15[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham"],:])
autocorplot(fit15[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham"],:])
out15 = DataFrame(fit15[501:end, :, :])
#remove burnin and save
CSV.write("Code/Results/Fit_BaseModel_nophoto.csv", out15)

#read in posterior predictiver
postpred15 = h5open("Code/Results/PostPred_BaseModel_nophoto.h5", "r") do f
    read(f, Chains)
end
postpredout15 = DataFrame(postpred15[501:end, :, :])

fit16 = h5open("Code/Results/Fit_BaseModel_justSOWER.h5", "r") do f
    read(f, Chains)
end
plot(fit16[501:end, ["r", "lnN1973", "K", "addCV_Branch"],:])
autocorplot(fit16[501:end, ["r", "lnN1973", "K", "addCV_Branch"],:])
out16 = DataFrame(fit16[501:end, :, :])

#remove burnin and save
CSV.write("Code/Results/Fit_BaseModel_justSOWER.csv", out16)

#read in posterior predictiver
postpred16 = h5open("Code/Results/PostPred_BaseModel_justSOWER.h5", "r") do f
    read(f, Chains)
end
postpredout16 = DataFrame(postpred16[501:end, :, :])
#save posterior predictive
CSV.write("Code/Results/PostPred_BaseModel_justSOWER.csv", postpredout16)



fit_11 = h5open("Code/Results/Fit_BaseModel_fixed_s.h5", "r") do f
    read(f, Chains)
end
plot(fit_11[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham"],:])
autocorplot(fit_11[501:end, ["r", "lnN1973", "K", "q_Ham", "addCV_Branch", "addCV_Ham"],:])
out_11 = DataFrame(fit_11[501:end, :, :])

#remove burnin and save
CSV.write("Code/Results/Fit_BaseModel_fixed_s.csv", out_11)
#read in posterior predictiver
postpred_11 = h5open("Code/Results/PostPred_BaseModel_fixed_s.h5", "r") do f
    read(f, Chains)
end
postpredout_11 = DataFrame(postpred_11[501:end, :, :])
#save posterior predictive
CSV.write("Code/Results/PostPred_BaseModel_fixed_s.csv", postpredout_11)
