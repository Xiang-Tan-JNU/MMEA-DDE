classdef MMEA_DDE < ALGORITHM
% <multi> <real/integer> <multimodal>
    methods
        function main(Algorithm, Problem)
            params.N = Problem.N;
            params.MaxIt = 50;
            params.L = 30;
            params.epsilon_b = 0.2;
            params.epsilon_n = 0.006;
            params.alpha = 0.5;
            params.delta = 0.995;
            params.T = 30;
            genFlag = [];

            netInitialized = 0;
            Pop1get = 0;

            eps = Algorithm.ParameterSet(0.3);
            flag_Ini_Arc=0;

            %% Generate random population
            Population          = Problem.Initialization();
            [~,CrowdDis1]       = EnvironmentalSelection(Population,Problem.N);
            [Archive,CrowdDis2,net,genFlag,flag_Ini_Arc] = ArchiveUpdate(Population,Problem.N,eps,0,Problem,params,[],genFlag,flag_Ini_Arc);

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                %% Initial the GNG network when can
                if ~netInitialized
                    NDNum = length(Archive);
                    if NDNum >= 2
                        net = InitilizeGrowingGasNet(Archive,params);
                        netInitialized = 1;
                    end
                end
                
                if Problem.FE <= Problem.maxFE*0.3 || ~netInitialized
                    MatingPool1 = TournamentSelection(2,round(Problem.N),-CrowdDis1);
                    Offspring   = OperatorGA(Problem,Population(MatingPool1));
                    [Population,CrowdDis1] = EnvironmentalSelection([Population,Offspring],Problem.N);
                    [Archive,CrowdDis2,net,genFlag,flag_Ini_Arc]    = ArchiveUpdate([Archive,Offspring],Problem.N,eps,Problem.FE/Problem.maxFE,Problem,params,net,genFlag,flag_Ini_Arc);
                else
                    if Pop1get == 0
                        Population1 = Archive;
                        Fitness1 = CalFitnessSup(Population1.decs,net.w);
                        Pop1get = 1;
                    end

                    MatingPool1 = TournamentSelection(2,Problem.N,-CrowdDis2);
                    Offspring1  = OperatorGAhalf(Problem,Archive(MatingPool1));

                    V = net.w;
                    MatingPool2 = randi(length(V),1,Problem.N);
                    Offspring2 = Problem.Evaluation(OperatorGAhalf(Problem,V(MatingPool2,:)));

                    MatingPool3 = TournamentSelection(2,round(Problem.N),-Fitness1);
                    Offspring3  = OperatorGAhalf(Problem,Population1(MatingPool3));

                    MatingPool4 = TournamentSelection(2,round(Problem.N),-CrowdDis1);
                    Offspring4  =  OperatorGAhalf(Problem,Population(MatingPool4));
                    Offspring=[Offspring1,Offspring2,Offspring3,Offspring4];

                    [Population,CrowdDis1] = EnvironmentalSelection([Population,Offspring],Problem.N);
                    [Population1,Fitness1] = EnvironmentalSelectionSup([Population1,Offspring],Problem.N,net);
                    [Archive,CrowdDis2,net,genFlag,flag_Ini_Arc]    = ArchiveUpdate([Archive,Offspring],Problem.N,eps,Problem.FE/Problem.maxFE,Problem,params,net,genFlag,flag_Ini_Arc);
                end
            end
        end
    end
end