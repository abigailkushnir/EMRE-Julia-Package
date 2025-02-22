using Pkg
Pkg.activate(".")
using LinearNoiseApproximation # need to import package from Julia to have LNA() function work
using DifferentialEquations
using Catalyst
using Reexport
using ModelingToolkit
using DiffEqBase
using JumpProcesses
using Plots
include("utils.jl") # need to import this utility file locally for LNA helper functions to work

using DifferentialEquations.EnsembleAnalysis
using Symbolics
using LinearAlgebra
using RowEchelon
using Statistics



## go through each reaction and adjust parameters to make units 1/t, i.e. multiply 0 --> X reactions by volume and divide biomolecular ones by volume, check if biomolecular reaction is same-species and if so multiply by 2
function adjustParams(rn::ReactionSystem, params, volume)

    systemParams = rn.ps
    systemReactions = rn.rxs
    
    adjustedTF = falses(length(systemParams))

    adjustedParams = copy(Float64.(params))

    for reaction in systemReactions
        if length(reaction.substrates) == 0 ## check if it's a zero order reaction
            paramIndex = findfirst(x->isequal(x,reaction.rate), systemParams) ## if it is, find what number parameter it is
            if !adjustedTF[paramIndex] ## if it hasn't already been adjusted
                adjustedParams[paramIndex] = adjustedParams[paramIndex]*volume ## multiply by volume
                adjustedTF[paramIndex] = true ## mark that it was adjusted
            end
        end
        if length(reaction.substrates) == 2 ## check if it's a biomolecular reaction
            paramIndex = findfirst(x->isequal(x,reaction.rate), systemParams) ## if it is, find what number parameter it is
            if !adjustedTF[paramIndex] ## if it hasn't already been adjusted
                adjustedParams[paramIndex] = adjustedParams[paramIndex]/volume ## divide by volume
                adjustedTF[paramIndex] = true
            end
        end
        if (length(reaction.substrates) == 1) && (reaction.substoich[1] == 2) ## check if it's a biomolecular reaction with same-species
            paramIndex = findfirst(x->isequal(x,reaction.rate), systemParams) ## if it is, find what number parameter it is
            if !adjustedTF[paramIndex] ## if it hasn't already been adjusted
                adjustedParams[paramIndex] = adjustedParams[paramIndex] * (2/volume) ## if it is, multiplay rate by 2/volume
                adjustedTF[paramIndex] = true ## mark that it's been adjusted
            end
        end
    end

    return adjustedParams
end;




function setParams(rn::ReactionSystem, params, volume)

    LNAparams = params;
    SSAparams = adjustParams(rn, params, volume) ## needs to know which rates are involved in biomolecular reactions as well as which ones need to be adjusted for units

    return LNAparams, SSAparams
end;






function runSSA(rn::ReactionSystem, ics, params, volume, tEnd, numSteps, numTrajs)

    tSpan = (0.0, tEnd); ## make time span
    timePoints = LinRange(0, tEnd, numSteps); ## set up time points to calculate SSA means at
    
    ## set up and run ensemble problem
    jumpsys = convert(JumpSystem, rn);
    dprob = DiscreteProblem(jumpsys, ics, tSpan, params);
    jprob = JumpProblem(jumpsys, dprob, Direct(), save_positions=(false, false));
    eprob = EnsembleProblem(jprob) 

    jsol = solve(eprob, SSAStepper(), trajectories = numTrajs, saveat = tEnd/(numSteps-1));

    means = timepoint_mean(jsol, timePoints)/volume; ## caclulate mean concentrations through time

    return timePoints, means
end;





function runLNA(rn::ReactionSystem, ics, params, tEnd)

    tSpan = (0.0, tEnd); ## set up time span

    ## set up and solve LNA problem from LinearNoiseApproximation.jl package
    lna = LNASystem(rn);
    prob = ODEProblem(lna, ics, tSpan, params);
    sol = solve(prob, Vern7(), saveat = tEnd/(numSteps-1));

    ## get indices of means and variances for each species in system
    meanIdxs, varIdxs = find_states_cov_number(collect(1:length(rn.species)), lna);

    return meanIdxs, varIdxs, sol

end;







function buildDeltaVect(rn::ReactionSystem, params, meanIdxs, varIdxs, LNAsol)

    deltaVec = zeros(length(rn.species))
    
    for reaction in rn.rxs

        if length(reaction.substrates) == 2 ## if a biomolecular reaction of different species
            species1idx = findfirst(x->isequal(x,reaction.substrates[1]), rn.species) ## find first species
            species2idx = findfirst(x->isequal(x,reaction.substrates[2]), rn.species) ## find second species

            covIdx = varIdxs[min(species1idx, species2idx)] + (max(species1idx, species2idx) - min(species1idx, species2idx)) ## find index of covariance, e.g. species 5 & 3 gives varIdxs[3] (which is variance of s3) + 2
            cov = LNAsol[covIdx, end] ## get covariance in SS (can change later to be in time)
            
            rateParamIdx = findfirst(x->isequal(x, reaction.rate), rn.ps) ## find what parameter is the rate 
            rate = params[rateParamIdx] ## get numerical rate from parameters
            
            adjustment = rate*cov ## k*<e_i,e_j>

            ## subtract the adjustment term from the two reactants
            deltaVec[species1idx] -= adjustment
            deltaVec[species2idx] -= adjustment
            
            for species in reaction.products ## for each product made in reaction

                prodIdx = findfirst(x->isequal(x,species),reaction.products) ## see what index it is in product list
                numProduced = reaction.prodstoich[prodIdx] ## see how many are produced by reaction

                speciesIdx = findfirst(x->isequal(x, species), rn.species) ## find index in the system species
                deltaVec[speciesIdx] += numProduced * adjustment ## add adjustment term for each one produced (??)
            end
        end

        if length(reaction.substrates) == 1 && (reaction.substoich[1] == 2) ## if it's a biomolecular reaction with the same species
            
            speciesIdx = findfirst(x->isequal(x, reaction.substrates[1]), rn.species) ## find what species it is

            ## find indices of its mean and variance
            meanIdx = meanIdxs[speciesIdx]
            varIdx = varIdxs[speciesIdx] 
            mean = LNAsol[meanIdx, end] ## SS mean (can change later to be in time)
            var = LNAsol[varIdx, end] ## SS var (can change later to be in time)

            rateParamIdx = findfirst(x->isequal(x, reaction.rate), rn.ps) ## find what parameter is the rate 
            rate = params[rateParamIdx] ## get numerical rate from parameters

            adjustment = rate * (var - mean) ## k * [<e_i^2> - phi_i]

            deltaVec[speciesIdx] -= 2 * adjustment ## subtract 2 * the adjustment from the reactants

            for species in reaction.products

                prodIdx = findfirst(x->isequal(x,species),reaction.products) ## see what index it is in product list
                numProduced = reaction.prodstoich[prodIdx] ## see how many are produced by reaction

                speciesIdx = findfirst(x->isequal(x, species), rn.species) ## find index in the system species
                deltaVec[speciesIdx] += numProduced * adjustment ## add adjustment term for each one produced (??)
            end
        end
    end

    return deltaVec
end;







function getInvJac(rn::ReactionSystem, params, meanIdxs, LNAsol)

    sys = convert(ODESystem, rn, combinatoric_ratelaws=false); ## convert reactioin system to ODE system
    jacMat = calculate_jacobian(sys); ## get the Jacobian
    
    subDict = Dict([]) ## initialize dictionary for substituting variables and steady states
    for species in rn.species ## for each species in the system
        speciesIdx = findfirst(x->isequal(x, species), rn.species) ## get index
        subDict[species] = LNAsol[meanIdxs[speciesIdx],end] ## if species shows up in Jacobian, substitute its steady state value
    end
    for param in rn.ps ## for each parameter in the system
        paramIdx = findfirst(x->isequal(x,param), rn.ps) ## get index
        subDict[param] = params[paramIdx] ## if parameter shows up in Jacobian, substitute the numerically value supplied
    end
    
    jacMatNumerical = substitute(jacMat, subDict) ## sub in numerical values into Jacobian
    jacMatInv = inv(jacMatNumerical) ## get inverse of Jacobian
    
    return jacMatInv
end;







function getRealizedVector(rn::ReactionSystem, params, meanIdxs, varIdxs, LNAsol, volume)

    SS = zeros(length(rn.species)) ## initialize vector for deterministic steady states
    for species in rn.species ## for each species
        speciesIdx = findfirst(x->isequal(x, species), rn.species) ## find index in system
        SS[speciesIdx] = LNAsol[meanIdxs[speciesIdx], end] ## record steady state value
    end

    invJac = getInvJac(rn, params, meanIdxs, LNAsol) ## get inverse Jacobian
    deltaVec = buildDeltaVect(rn, params, meanIdxs, varIdxs, LNAsol) ## get delta vector with adjustments

    realizedSS = SS - (1/volume)*invJac*deltaVec ## Eq. 26 from 2009 paper

    return SS, realizedSS
end;







function EMRE(rn::ReactionSystem, ics, params, volume, tEnd, numSteps, numTrajs)

    LNAparams, SSAparams = setParams(rn, params, volume)
 
    timePoints, SSAmeans = runSSA(rn, ics, SSAparams, volume, tEnd, numSteps, numTrajs)
    meanIdxs, varIdxs, LNAsol = runLNA(rn, ics, LNAparams, tEnd)
 
    SS, realizedSS = getRealizedVector(rn, LNAparams, meanIdxs, varIdxs, LNAsol, volume)
    
    return timePoints, SSAmeans, SS, realizedSS
 end;



