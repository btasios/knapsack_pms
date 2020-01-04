#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include "../ortools/include/ortools/algorithms/knapsack_solver.h"
#include "../ortools/include/ortools/linear_solver/linear_solver.h"


using namespace std;

struct item
{
    int id;
    int profit;
    int weight;
};

struct knapsack_problem
{
    vector<item> items;
    int capacity;
    int time_limit = 10; // seconds
};

struct bnbNode{
    vector<item> items_included;
    vector<item> items_remaining;
    int profit;
    int weight;
    int upper_bound;
};

struct solution{
    vector<item> items;
    double time;
};

knapsack_problem read_data(string &fn)
{
    ifstream fin(fn);
    if (!fin.good())
    {
        cerr << "Error opening file " << fn << endl;
        exit(-1);
    }
    knapsack_problem ks;
    int items_number;
    fin >> items_number;
    for (int i = 0; i < items_number; i++)
    {
        item an_item;
        fin >> an_item.id;
        fin >> an_item.profit;
        fin >> an_item.weight;
        ks.items.push_back(an_item);
    }
    fin >> ks.capacity;
    return ks;
}

void print_knapsack_problem_info(knapsack_problem &ks)
{
    cout << "Items=" << ks.items.size() << std::endl;
    for (int i = 0; i < ks.items.size(); i++)
    {
        std::cout << "Item=" << ks.items[i].id << " value=" << ks.items[i].profit << " weight=" << ks.items[i].weight << std::endl;
    }
    std::cout << "Capacity=" << ks.capacity << std::endl;
}

int get_profit(knapsack_problem &ks, vector<item> &sol){ // Will be used Later!
    int total_weight = 0;
    int total_profit = 0;
    for(item an_item:sol){
        total_weight += an_item.weight;
        if(total_weight > ks.capacity) return -1;
        total_profit += an_item.profit;
    }
    return total_profit;
}

void export_solution(knapsack_problem &ks, vector<item> &sol, string fn){
    ofstream ofs(fn);
    if (!ofs.good()){
        cerr << "An error accessing the file has occured!" << endl;
        exit(-1);
    }
    string output;
    int total_weight=0;
    int total_profit=0;
    for(int i=0; i<sol.size();i++){
        output=output + to_string(sol[i].id) + " ";
        total_weight += sol[i].weight;
        total_profit += sol[i].profit;
    }
    output += to_string(total_profit) + " " + to_string(total_weight);
    ofs << output;
    ofs.close();
}

void export_to_csv(knapsack_problem &ks, solution &a_sol, ofstream &mycsv){
    string output;
    int total_weight=0;
    int total_profit=0;
    for(int i=0; i<a_sol.items.size();i++){
        //output=output + to_string(a_sol.items[i].id) + "-";
        total_weight += a_sol.items[i].weight;
        total_profit += a_sol.items[i].profit;
    }
    output =to_string(total_profit) + "," + to_string(total_weight) + "," + to_string(a_sol.time) + ",";
    mycsv << output;

}

void add_headers(ofstream &mycsv){
    mycsv << "Problem" << "," << "Greedy Profit" << "," << "Greedy Weight" << "," << "Greedy Time" << ",";
    mycsv << "Brute Force Profit" << "," << "Brute Force Weight" << "," << "Brute Force Time" << ",";
    mycsv << "Branch and Bound Profit" << "," << "Branch and Bound Weight" << "," << "Branch and Bound Time" << ",";
    mycsv << "Dynamic Profit" << "," << "Dynamic Weight" << "," << "Dynamic Time" << ",";
    mycsv << "OR_Dynamic Profit" << "," << "OR_Dynamic Weight" << "," << "OR_Dynamic Time" << ",";
    mycsv << "OR_Integer Profit" << "," << "OR_Integer Weight" << "," << "OR_Integer Time" << ",";
    mycsv << "\n";
}


solution greedy_solver(knapsack_problem &ks){ // H lysh me Greedy
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    vector<item> itemsQueue = ks.items;
    sort(itemsQueue.begin(),itemsQueue.end(), [](item &item1, item &item2){ return (double) item1.profit/(double)item1.weight>(double)item2.profit/(double)item2.weight;});

    int total_weight =0;
    solution a_solution;
    for(int i=0; i<itemsQueue.size();i++){
        if(total_weight+itemsQueue[i].weight>ks.capacity) break;
        a_solution.items.push_back(itemsQueue[i]);
        total_weight += itemsQueue[i].weight;
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
        if(time_span.count()>ks.time_limit){
            cout << "Greedy: Time limit has been reached!" << endl;
            goto TIMED;
        }
    }
    TIMED:chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    sort(a_solution.items.begin(),a_solution.items.end(), [](item &item1, item &item2){return item1.id < item2.id;});
    cout << "Time for Greedy:" << std::fixed << time_span.count() << " seconds" << endl;
    a_solution.time = time_span.count();
    return a_solution;
}

solution brute_force(knapsack_problem &ks){ // H lysh me Brute force
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    vector<vector<item>> powerset;
    solution a_solution;
    int max_profit = -1;
    for(item an_item: ks.items){
        vector<vector<item>> new_sets;
        for(vector<item> a_set : powerset){
            a_set.push_back(an_item);
            new_sets.push_back(a_set);
            int thisprofit = get_profit(ks,a_set);
            if(thisprofit>max_profit){
                max_profit=thisprofit;
                a_solution.items=a_set;
            }
            chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
            chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
            if(time_span.count()>ks.time_limit){
                cout << "Brute Force: Time limit has been reached!" << endl;
                cout << "This time is:" << std::fixed << time_span.count() << " seconds" << endl;
                goto TIMED;
            }
        }
        
        powerset.insert(powerset.end(), new_sets.begin(), new_sets.end());
        powerset.push_back({an_item});
    }
    TIMED:powerset.push_back({}); // bazoume to keno
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cout << "Time for Brute Force:" << std::fixed << time_span.count() << " seconds" << endl;
    sort(a_solution.items.begin(),a_solution.items.end(), [](item &item1, item &item2){return item1.id < item2.id;});
    a_solution.time=time_span.count();
    return a_solution;   
}

int bound(knapsack_problem &ks, vector<item> &items_remaining,int &current_profit, int &current_weight){
    int weight_remaining=ks.capacity-current_weight;
    int upper_bound = current_profit;
    int i = 0;
    while(weight_remaining>0){
        if(weight_remaining-items_remaining[i].weight>=0){
            weight_remaining=weight_remaining-items_remaining[i].weight;
            upper_bound=upper_bound+items_remaining[i].profit;
        }else{
            upper_bound=upper_bound+int(items_remaining[i].profit/items_remaining[i].weight)*weight_remaining;
            weight_remaining=0;
        }
        i++;
    }
    return upper_bound;
}

void find_current_max(solution &a_solution, vector<bnbNode> &liveNodes,int &max_profit){
    for(bnbNode node: liveNodes){
        if(max_profit<node.profit){
            max_profit=node.profit;
            a_solution.items=node.items_included;
        } 
    }
}

solution branch_and_bound(knapsack_problem &ks){
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    vector<item> itemsQueue = ks.items;
    sort(itemsQueue.begin(),itemsQueue.end(), [](item &item1, item &item2){ return (double) item1.profit/(double)item1.weight>(double)item2.profit/(double)item2.weight;});
    bnbNode root,left,right;
    vector<bnbNode> liveNodes;
    solution a_solution;
    root.profit=0;
    root.weight=0;
    root.items_remaining=itemsQueue;
    root.upper_bound=bound(ks,root.items_remaining,root.profit,root.weight);
    liveNodes.push_back(root);
    int max_profit=0;
    while (!liveNodes.empty())
    {
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
        if(time_span.count()>ks.time_limit){
            cout << "Branch and bound: Time limit has been reached!" << endl;
            goto TIMED;
        }
        sort(liveNodes.begin(),liveNodes.end(), [](bnbNode &n1, bnbNode &n2){return (double) n1.upper_bound>n2.upper_bound;}); // sort by upper bound
        find_current_max(a_solution,liveNodes,max_profit);
        /*cout << "Solution:" << endl;
        for(item an_item: a_solution.items){
            cout << an_item.id << " ";
        }
        cout << endl;*/
        if(!liveNodes.front().items_remaining.empty() && liveNodes.front().upper_bound>max_profit && liveNodes.front().profit>=0){
            left.items_included=liveNodes.front().items_included;
            left.items_included.push_back(liveNodes.front().items_remaining[0]);
            right.items_included=liveNodes.front().items_included;
            left.profit=get_profit(ks,left.items_included);
            right.profit=get_profit(ks,right.items_included);
            left.weight=liveNodes.front().weight+liveNodes.front().items_remaining[0].weight;
            right.weight=liveNodes.front().weight;
            liveNodes.front().items_remaining.erase(liveNodes.front().items_remaining.begin());
            left.items_remaining = liveNodes.front().items_remaining;
            //left.items_remaining.erase(left.items_remaining.begin());
            right.items_remaining= liveNodes.front().items_remaining;
            //right.items_remaining.erase(right.items_remaining.begin());
            left.upper_bound=bound(ks,left.items_remaining,left.profit,left.weight);
            right.upper_bound=bound(ks,right.items_remaining,right.profit,right.weight);
            liveNodes.erase(liveNodes.begin());
            liveNodes.push_back(left);
            liveNodes.push_back(right);
        }else{
            liveNodes.erase(liveNodes.begin());
        }
    }
    TIMED:chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cout << "Time for Branch and Bound:" << std::fixed << time_span.count() << " seconds" << endl;
    sort(a_solution.items.begin(),a_solution.items.end(), [](item &item1, item &item2){return item1.id < item2.id;});
    a_solution.time=time_span.count();
    //Sort it, time it, change it to solution struct
    return a_solution;
}

solution dynamic(knapsack_problem &ks){ // H lysh me Dynamic
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    int n=ks.items.size();
    int C=ks.capacity;
    int** V = new int* [n+1];
    for(int i=0;i<=n;i++)
        V[i] = new int[C+1];
    //int V[n+1][C+1];
    for(int j=0;j<=C;j++)
        V[0][j] = 0;
    for(int i=1; i<=n;i++){ // opou yparxei ks.items einai -1 epeidh to item id:1 einai sth thesh 0 tou ks.items kok
        for(int j=0;j<=C;j++){
            if(ks.items[i-1].weight > j){
                V[i][j] = V[i-1][j];
            }else{
                V[i][j] = max(V[i-1][j], V[i-1][j-ks.items[i-1].weight] + ks.items[i-1].profit); 
            }
            chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
            chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
            if(time_span.count()>ks.time_limit){
                cout << "Dynamic: Time limit has been reached!" << endl;
                goto TIMED;
            }
        }
    }
    TIMED:solution a_solution;
    int j = C;
    for(int i=n;i>0;i--)
        if((ks.items[i-1].weight<=j) && ((V[i-1][j-ks.items[i-1].weight]+ks.items[i-1].profit) >= V[i-1][j])){
            a_solution.items.push_back(ks.items[i-1]);
            j = j - ks.items[i-1].weight;
        }
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cout << "Dynamic time:" << std::fixed << time_span.count() << " seconds" << endl;
    for(int i=0;i<=n;i++)
        delete[] V[i];
    delete[] V;
    sort(a_solution.items.begin(),a_solution.items.end(), [](item &item1, item &item2){return item1.id < item2.id;});
    a_solution.time = time_span.count();
    return a_solution;
}

namespace operations_research
{
solution dynamic_ortools(knapsack_problem &ks){
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    KnapsackSolver solver(KnapsackSolver::KNAPSACK_DYNAMIC_PROGRAMMING_SOLVER,
                          "SimpleKnapsackExample");
    std::vector<int64> values;
    std::vector<std::vector<int64>> weights;
    std::vector<int64> temp_weights;
    std::vector<int64> capacities = {ks.capacity};
    solution a_solution;
    for(item an_item: ks.items){
        values.push_back(an_item.profit);
        temp_weights.push_back(an_item.weight);
    }
    weights.push_back(temp_weights);

    solver.Init(values, weights, capacities);
    int64 computed_value = solver.Solve();
    
    std::vector<int> packed_items;
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (solver.BestSolutionContains(i))
            a_solution.items.push_back(ks.items[i]);        
    }
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    std::cout << "Time for ORTOOLS Dynamic: " << std::fixed << time_span.count() << " seconds" << endl;
    a_solution.time = time_span.count();
    return a_solution;
}

solution integer_ortools(knapsack_problem &ks){
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    MPSolver solver("IP KNAPSACK SOLVER", MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);
    MPConstraint *const c = solver.MakeRowConstraint(0, ks.capacity, "capacity_constraint");
    MPObjective *const objective = solver.MutableObjective();
    std::vector<MPVariable*> items_x;                                                  
    for (int i=0;i<ks.items.size();i++){
        std::string s= "x" + std:: to_string((i+1));
        //std::cout<<s<<"\n"; 
        MPVariable *const x = solver.MakeBoolVar(s);
        items_x.push_back(x);
        c->SetCoefficient(x, ks.items.at(i).weight);
        objective->SetCoefficient(x, ks.items.at(i).profit);
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
            if(time_span.count()>ks.time_limit){
                cout << "Integer: Time limit has been reached!" << endl;
                goto TIMED;
            }
    }
    TIMED:objective->SetMaximization();
    solver.Solve();
    solution a_solution;
    for (int i=0;i<items_x.size();i++){
       int in_solution=items_x.at(i)->solution_value();
       if(in_solution==1) a_solution.items.push_back(ks.items.at(i));
    }
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    std::cout << "Time for ORTOOLS Integer: " << std::fixed << time_span.count() << " seconds" << endl;
    a_solution.time = time_span.count();
    return a_solution;
}
} // namespace operations_research

int main()
{
    ofstream mycsv;
    if (!mycsv.good()){
        cerr << "An error accessing the csv has occured!" << endl;
        exit(-1);
    }
    mycsv.open("results.csv");
    add_headers(mycsv);
    for(int i:{1, 2, 3, 4, 5})
        for(int n: {10, 50, 100, 500})
            for(int r: {50, 100, 500, 1000})
                for(int t: {1, 2, 3, 4}){
                    string fn = "problem_" + to_string(n) + "_"+ to_string(r) + "_" + to_string(t) + "_" + to_string(i) + "_5.txt";
                    knapsack_problem ks = read_data(fn);
                    mycsv << fn + ",";
                    solution sol;
                    sol = greedy_solver(ks);
                    export_to_csv(ks,sol,mycsv);
                    export_solution(ks,sol.items,"GR_"+fn);
                    sol = brute_force(ks);
                    export_to_csv(ks,sol,mycsv);
                    export_solution(ks,sol.items,"BF_"+fn);
                    sol = branch_and_bound(ks);
                    export_to_csv(ks,sol,mycsv);
                    export_solution(ks,sol.items,"BNB_"+fn);
                    sol = dynamic(ks);
                    export_to_csv(ks,sol,mycsv);
                    export_solution(ks,sol.items,"D_"+fn);
                    sol=operations_research::dynamic_ortools(ks);
                    export_to_csv(ks,sol,mycsv);
                    export_solution(ks,sol.items,"OR_D_"+fn);
                    sol=operations_research::integer_ortools(ks);
                    export_to_csv(ks,sol,mycsv);
                    export_solution(ks,sol.items,"OR_I_"+fn);

                    mycsv << "\n";
                }

    mycsv.close();
}