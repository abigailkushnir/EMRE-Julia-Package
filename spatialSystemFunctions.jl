include("EMREfunctionsCompleteSS.jl");



function createSubSpecies(symbol::Symbol, voxelNum::Int)
    subSymbol =  string(string(symbol)[1], string(voxelNum))
    return Symbol(subSymbol)
end;





function createSubRate(symbol::Symbol, voxelNum::Int)
    subSymbol =  string(string(symbol), string(voxelNum))
    return Symbol(subSymbol)
end;





function createSubDiffRate(symbol::Symbol)
    subSymbol =  string("d", string(symbol))
    return Symbol(subSymbol)
end;





struct VoxelReaction
    rxn::Reaction
    voxelNums::Vector{Int}
    voxelRates::Vector{Float64}
end





struct VoxelSpecies
    speciesName::Num
    voxelsNums::Vector{Int}
    initConds::Vector{Float64}
end;





struct DiffReaction
    speciesName::Num
    voxelPairs::Vector{}
    voxelRates::Vector{Float64}
end;





@kwdef struct SpatialSystem
    sysSpecies::Vector{VoxelSpecies}
    sysReactions::Vector{VoxelReaction}
    sysDiffReactions::Vector{DiffReaction}
    numVoxels::Int
end;





function addSpecies(rn::ReactionSystem, species::VoxelSpecies, icsDict)
    i = 1;
    for voxelNum in species.voxelsNums
        subSpeciesSymb = createSubSpecies(Symbol(species.speciesName), voxelNum)
        subSpecies = (@species ($subSpeciesSymb)(t))[1]
        
        addspecies!(rn, subSpecies)
        icsDict[subSpecies] = species.initConds[i]
        i += 1
    end;
end;





function addSystemReaction(rn::ReactionSystem, reaction::VoxelReaction, paramsDict)
    i = 1;
    for voxelNum in reaction.voxelNums
        voxelSubstrates = [];
        for substrate in reaction.rxn.substrates
            subSpeciesSymb = createSubSpecies(Symbol(substrate), voxelNum)
            subSpecies = (@species ($subSpeciesSymb)(t))[1]
            push!(voxelSubstrates, subSpecies)
            addspecies!(rn, subSpecies)
        end;
        voxelProducts = [];
        for product in reaction.rxn.products
            prodSpeciesSymb = createSubSpecies(Symbol(product), voxelNum)
            prodSpecies = (@species ($prodSpeciesSymb)(t))[1]
            push!(voxelProducts, prodSpecies)
            addspecies!(rn, prodSpecies)
        end;

        if length(reaction.voxelRates) == 1
            react = Reaction(reaction.rxn.rate, voxelSubstrates, voxelProducts, reaction.rxn.substoich, reaction.rxn.prodstoich)
            addreaction!(rn, react)
            addparam!(rn, reaction.rxn.rate)
            paramsDict[reaction.rxn.rate] = reaction.voxelRates[1]
        else
            rateSymb = createSubRate(Symbol(reaction.rxn.rate), voxelNum)
            rate = (@parameters $rateSymb)[1]
            react = Reaction(rate, voxelSubstrates, voxelProducts, reaction.rxn.substoich, reaction.rxn.prodstoich)
            addreaction!(rn, react)
            addparam!(rn, rate)
            paramsDict[rate] = reaction.voxelRates[i]
        end;
        i += 1
    end;
    addparam!(rn, reaction.rxn.rate); ## assuming rate is the same across voxels, move up into for loop if this changes to different rates
end;





## adds diffusion reactions for a species to appropriate voxels
function addSpeciesDiffusion(rn::ReactionSystem, diff::DiffReaction, paramsDict)
    i = 1;
    for pair in diff.voxelPairs
        subSpeciesSymb1 = createSubSpecies(Symbol(diff.speciesName), pair[1]);
        subSpecies1 = (@species ($subSpeciesSymb1)(t))[1];
        subSpeciesSymb2 = createSubSpecies(Symbol(diff.speciesName), pair[2]);
        subSpecies2 = (@species ($subSpeciesSymb2)(t))[1];
    
        addspecies!(rn, subSpecies1);
        addspecies!(rn, subSpecies2);

        rateSymb = createSubDiffRate(Symbol(diff.speciesName));

        if length(diff.voxelRates) == 1
            rate = (@parameters $rateSymb)[1]
            react1 = Reaction(rate, [subSpecies1], [subSpecies2]);
            addreaction!(rn, react1);
            addparam!(rn, rate);
            paramsDict[rate] = diff.voxelRates[1]

            if pair[3] == true ## if diffusion is bidirectional
                react2 = Reaction(rate, [subSpecies2], [subSpecies1]);
                addreaction!(rn, react2);
            end;
        else
            ## can't do this yet
        end;
    end;
end;





function buildSystem(sys::SpatialSystem)

    rn = @reaction_network begin
    end;

    icsDict = Dict([]);
    paramsDict = Dict([]);

    for species in sys.sysSpecies
        addSpecies(rn, species, icsDict)
    end;
    for rx in sys.sysReactions
        addSystemReaction(rn, rx, paramsDict)
    end;
    for dRx in sys.sysDiffReactions
        addSpeciesDiffusion(rn, dRx, paramsDict)
    end;

    return rn, icsDict, paramsDict;
end;





## a function that builds the spatial system and sets up the ics and params for EMRE functions
function getSysICsParams(sys::SpatialSystem)

    rn, icsDict, paramsDict = buildSystem(sys);
    
    ics = [];
    for species in rn.species
        push!(ics, icsDict[species]);
    end;
    params = [];
    for param in rn.ps
        if param in keys(paramsDict)
            push!(params, paramsDict[param]);
        else
            push!(params, 0)
        end;
    end;

    return rn, ics, params;
end;





