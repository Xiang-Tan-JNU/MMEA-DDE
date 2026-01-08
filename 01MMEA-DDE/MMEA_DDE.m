classdef MMEA_DDE < ALGORITHM
% <multi> <real/integer> <multimodal>
% eps --- 0.3 --- 
    methods
        function main(Algorithm, Problem)
            % Parameters for GNG network
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

            eps = Algorithm.ParameterSet(0.3);
            flag_Ini_Arc=0;

            Population          = Problem.Initialization();
            [~,CrowdDis1]       = EnvironmentalSelection(Population,Problem.N);
            [Archive,CrowdDis2,net,genFlag,flag_Ini_Arc] = ArchiveUpdate(Population,Problem.N,eps,0,Problem,params,[],genFlag,flag_Ini_Arc);

            while Algorithm.NotTerminated(Archive)
                if ~netInitialized
                    NDNum = length(Archive);
                    if NDNum >= 2
                        net = InitilizeGrowingGasNet(Archive,params);
                        netInitialized = 1;
                    end
                end
                if Problem.FE <= Problem.maxFE*0.2 || ~netInitialized
                    MatingPool1 = TournamentSelection(2,round(Problem.N),-CrowdDis1);
                    Offspring   = OperatorGA(Problem,Population(MatingPool1));
                    [Population,CrowdDis1] = EnvironmentalSelection([Population,Offspring],Problem.N);
                    [Archive,CrowdDis2,net,genFlag,flag_Ini_Arc]    = ArchiveUpdate([Archive,Offspring],Problem.N,eps,Problem.FE/Problem.maxFE,Problem,params,net,genFlag,flag_Ini_Arc);
                else
                    MatingPool1 = TournamentSelection(2,Problem.N,-CrowdDis2);
                    Offspring1  = OperatorGAhalf(Problem,Archive(MatingPool1));

                    V = net.w;
                    MatingPool2 = randi(length(V),1,Problem.N);
                    Offspring2 = Problem.Evaluation(OperatorGAhalf(Problem,V(MatingPool2,:)));

                    MatingPool4 = TournamentSelection(2,round(Problem.N),-CrowdDis1);
                    Offspring4  =  OperatorGAhalf(Problem,Population(MatingPool4));
                    Offspring=[Offspring1,Offspring2,Offspring4];

                    [Population,CrowdDis1] = EnvironmentalSelection([Population,Offspring],Problem.N);
                    [Archive,CrowdDis2,net,genFlag,flag_Ini_Arc]    = ArchiveUpdate([Archive,Offspring],Problem.N,eps,Problem.FE/Problem.maxFE,Problem,params,net,genFlag,flag_Ini_Arc);
                end
            end
        end
    end
end
