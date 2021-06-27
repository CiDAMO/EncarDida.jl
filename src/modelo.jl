using CSV: DataFrames
using JuMP, CSV, DataFrames, Gurobi, LinearAlgebra, Random, XLSX, StatsBase

StatsBase.rand 
Random.seed!(0)

#Dados
lista_materias = DataFrame(XLSX.readtable("MateriasFinais.xlsx", "materiasFinais")...)

lista_professores = DataFrame(XLSX.readtable("preferenciasFinais.xlsx", "preferenciasFinais")...)

lista_restricoes = DataFrame(XLSX.readtable("RestricoesFinais.xlsx", "RestricoesFinais")...)

function FacultyAssign( DataFrameMaterias ,
                        DataFramePreferencias,
                        DataFrameRestricoes;
                        ListaProfesPos = [6, 49, 4, 0, 38, 22, 19, 44, 45, 40, 0, 27, 0],
                        ListaMateriasPos = 93:105,
                        ListaProfesComCargos = [52, 20, 30, 46],
                        ListaProfesSubs = [42, 53, 54],
                        NumeroHorarios = 8,
                        NumeroDias = 5,
                        CargaHorariaComumMin = 8,
                        CargaHorariaComumMax = 12,
                        CargaHorariaCargoMin = 4,
                        CargaHorariaCargoMax = 6,
                        CargaHorariaSubsMin = 12,
                        CargaHorariaSubsMax = 16
                        )


    P = size(DataFramePreferencias, 1) # número de professores e matérias
    M = size(DataFramePreferencias, 2) - 1 # número de professores e matérias
    H = NumeroHorarios # horarios 
    D = NumeroDias # dias

    # Gerando matriz de duração de aulas 
    DU = zeros(D, H, M)
    for k = 1:M
        DU[:, :, k] .= DataFrameMaterias[1 + (k - 1) * 5: 5 + (k - 1) * 5, 2:(H + 1)]
    end

    # Gerando a matriz Binária de aula
    HT = copy(DU) / 2

    # Gerando a matriz das preferências
    preferencias = Matrix(DataFramePreferencias[:, 2:(M+1)])

    # Gerando a matriz de restrições 
    restricoes = zeros(D, H, P)
    for k = 1:P
        restricoes[:, :, k] .= DataFrameRestricoes[1 + (k - 1) * 5: 5 + (k - 1) * 5, 2: end]
    end


    ## Completando professores da pos 
    sorteioProfePos = sample(setdiff(1:P, ListaProfesPos), 3, replace = false)
    ListaProfesPos[4] = sorteioProfePos[1]
    ListaProfesPos[11] = sorteioProfePos[2]
    ListaProfesPos[13] = sorteioProfePos[3]

    preferencias[sorteioProfePos[1],96] = 5
    preferencias[sorteioProfePos[2],103] = 5
    preferencias[sorteioProfePos[3],105] = 5


    # Definindo o modelo
    model = Model(Gurobi.Optimizer)
    @variable(model, x[1:P, 1:M] ≥ 0, Bin) #P*M var de decisão

    @constraint(model,[p=1:P, d=1:D, h=1:H], sum(HT[d, h, t] * x[p,t] for t=1:M) ≤ 1)

    @constraint(model, [t=1:M], sum(x[p,t] for p=1:P) == 1)

    # Professores "comuns"
    @constraint(model, [p in setdiff(1:P, hcat(ListaProfesComCargos', ListaProfesSubs'))], sum(DU[d, h, t]*x[p,t] for t=1:M, d=1:D, h=1:H) ≥ CargaHorariaComumMin)
    @constraint(model, [p in setdiff(1:P, hcat(ListaProfesComCargos', ListaProfesSubs'))], sum(DU[d, h, t]*x[p,t] for t=1:M, d=1:D, h=1:H) ≤ CargaHorariaComumMax)

    # Professores com cargos 
    @constraint(model, [p in ListaProfesComCargos'], sum(DU[d, h, t]*x[p,t] for t=1:M, d=1:D, h=1:H) ≥ CargaHorariaCargoMin)
    @constraint(model, [p in ListaProfesComCargos'], sum(DU[d, h, t]*x[p,t] for t=1:M, d=1:D, h=1:H) ≤ CargaHorariaCargoMax)

    # Professores da Pós
    @constraint(model, [i=1:length(ListaMateriasPos)], x[ListaProfesPos[i], ListaMateriasPos[i]] == 1)

    # Professores Substitutos
    @constraint(model, [p in ListaProfesSubs], sum(DU[d, h, t]*x[p,t] for t=1:M, d=1:D, h=1:H) ≥ CargaHorariaSubsMin)
    @constraint(model, [p in ListaProfesSubs], sum(DU[d, h, t]*x[p,t] for t=1:M, d=1:D, h=1:H) ≤ CargaHorariaSubsMax)

    @constraint(model, [p=1:P, t=1:M], x[p,t] * sum(HT[d, h, t] * restricoes[d, h, p] for h=1:H, d=1:D) == 0)

    for t=setdiff(1:M, [24,85,86]) 
        if HT[2,5,t] == 1
            @constraint(model, [p=1:P], x[p, 85] + x[p, t] ≤ 1)
        end
        if HT[1,5,t] == 1
            @constraint(model, [p=1:P], x[p, 86] + x[p, t] ≤ 1)
        end
        if HT[2,6,t] == 1
            @constraint(model, [p=1:P], x[p, 24] + x[p, t] ≤ 1)
        end
    end

    @objective(model, Max, sum(preferencias[p,t]*x[p,t] for p=1:P, t=1:M));

    optimize!(model) #resolver


    solucao = value.(x)
    df = DataFrame([DataFramePreferencias[:,1] solucao])
    rename!(df, names(DataFramePreferencias))   


    SomaPesoProf = sum(solucao .* preferencias, dims=2)
    mediaSomaPref = sum(SomaPesoProf) / size(DataFramePreferencias, 1)
    desviosSomaPref = SomaPesoProf - ones(size(DataFramePreferencias, 1)) * mediaSomaPref
    sdPref = sqrt((transpose(desviosSomaPref) * desviosSomaPref) / length(desviosSomaPref))
    maxPref = maximum(SomaPesoProf)
    minPref = minimum(SomaPesoProf)
    Resultado1 = DataFrame([mediaSomaPref sdPref maxPref minPref],["média", "DesvioPadrao" ,"max", "min"]) 

    A = solucao .* preferencias
    MateriaPeso5 = length(findall(A .== 5))
    MateriaPeso4 = length(findall(A .== 4))
    MateriaPeso3 = length(findall(A .== 3))
    MateriaPeso2 = length(findall(A .== 2))
    MateriaPeso1 = length(findall(A .== 1))
    MateriaPeso0 = 92 - sum([MateriaPeso1 MateriaPeso2 MateriaPeso3 MateriaPeso4 MateriaPeso5])
    Resultado2 = DataFrame([MateriaPeso0 MateriaPeso1 MateriaPeso2 MateriaPeso3 MateriaPeso4 MateriaPeso5],["peso 0", "peso 1", "peso 2", "peso 3", "peso 4", "peso 5"]) 

    soma = sum([DataFramePreferencias[:,1] solucao][:,2:end], dims=2)
    materias3 = length(findall(soma .== 3.0)) # professores com 3 matérias
    materias2 = length(findall(soma .== 2.0)) # professores com 2 matérias
    materias1 = length(findall(soma .== 1.0)) # professores com 1 matéria
    quantidade_materias = DataFrame([materias3 materias2 materias1], ["3 matérias", "2 matérias", "1 matéria"])

    Status = hcat([Resultado1 Resultado2 quantidade_materias])

    return Dict(:DataFrameSolucao => df,
            :StatusGeral => Status,
            :Estatisticas => Resultado1,
            :DistribuicaoPesos => Resultado2,
            :QuantidadeMaterias => quantidade_materias
        )

end

resolve = FacultyAssign( lista_materias,
               lista_professores,
               lista_restricoes)


resolve[:StatusGeral]

resolve[:Estatisticas]

resolve[:DistribuicaoPesos]

resolve[:QuantidadeMaterias]