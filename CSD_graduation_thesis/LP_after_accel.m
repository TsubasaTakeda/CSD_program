per_driver_array = [64];
num_node_array = [49];
num_depot_array = [8];

num_networks = 1;
num_experience = 10;


for num_node = num_node_array
    for num_depot = num_depot_array
        for per_driver = per_driver_array
            
            sum_driver = per_driver*num_node*num_node;
            sum_task = sum_driver;
    
            
            
            filename = strcat("result/result_proposed/num_driver_", num2str(sum_driver), "/node_", num2str(num_node), "/depot_", num2str(num_depot));
            max_time_filename = strcat(filename, "/time/LP_max_time_", num2str(sum_driver), "_", num2str(num_node), "_", num2str(num_depot), ".csv");
            sum_time_filename = strcat(filename, "/time/LP_sum_time_", num2str(sum_driver), "_", num2str(num_node), "_", num2str(num_depot), ".csv");
            obj_filename = strcat(filename, "/obj/proposed_obj_", num2str(sum_driver), "_", num2str(num_node), "_", num2str(num_depot), ".csv");
            
            foldername = strcat(filename, "/obj");
            mkdir(foldername);
            
            
            proposal_obj_array = zeros(num_networks, num_experience);
            max_LP_time_array = zeros(num_networks, num_experience);
            sum_LP_time_array = zeros(num_networks, num_experience);
            
            for net_num = 0:num_networks-1
                for exp_num = 0:num_experience-1
                    
                    
                    filename_input = strcat("input/num_driver_", num2str(sum_driver), "/node_", num2str(num_node), "/depot_", num2str(num_depot) + "/net_", num2str(net_num));
                    filename_output = strcat("result/result_proposed/num_driver_", num2str(sum_driver), "/node_", num2str(num_node), "/depot_", num2str(num_depot), "/net_", num2str(net_num));
                    
                    %ODペアごとのドライバ�?�数?��タスク数?���?�別知覚費用を取得するため�?�リンク
                    driver_cost_atomic_filename = strcat(filename_input, "/cost_atomic/cost_atomic_driver_", num2str(num_node), "_", num2str(num_depot), "_", num2str(net_num), "_", num2str(exp_num), ".csv");
                    task_cost_atomic_filename = strcat(filename_input, "/cost_atomic/cost_atomic_task_", num2str(num_node), "_", num2str(num_depot), "_", num2str(net_num), "_", num2str(exp_num), ".csv");
                    demand_filename = strcat(filename_input, "/demand/demand_", num2str(num_node), "_", num2str(num_depot), "_", num2str(net_num), "_" + num2str(exp_num), ".csv");
                    
                    %�?速勾配法によって決定した�?��?を取得するため�?�リンク
                    driver_allocation_filename = strcat(filename_output, "/allocation/allocation_", num2str(exp_num), "/rounded_driver_allocation_", num2str(num_node), "_", num2str(num_depot), "_", num2str(net_num), "_", num2str(exp_num), ".csv");
                    task_allocation_filename = strcat(filename_output, "/allocation/allocation_", num2str(exp_num), "/rounded_task_allocation_", num2str(num_node), "_", num2str(num_depot), "_", num2str(net_num), "_", num2str(exp_num), ".csv");
                    
                    %コスト行�?�，ドライバ�?�数?��タスク数の取�?
                    driver_cost = readmatrix(driver_cost_atomic_filename);
                    task_cost = readmatrix(task_cost_atomic_filename);
                    demand = round(readmatrix(demand_filename));
                    
                    
                    %配�?を取�?
                    driver_allocation = readmatrix(driver_allocation_filename);
                    task_allocation = readmatrix(task_allocation_filename);
                    
                    %ODペアごとのドライバ�?�数?��タスク数
                    num_driver = demand(1:num_node^2).';
                    num_task = demand(num_node^2+1:num_node^2 + num_depot*num_node).';
                    
                    clear filename_input;
                    clear filename_output;
                    clear driver_cost_atomic_filename;
                    clear task_cost_atomic_filename;
                    clear driver_allocation_filename;
                    clear task_allocation_filename;
                    clear demand_filename;
                    clear demand;
                    
                    %driver_cost = reshape(driver_cost.', 1, []);
                    %task_cost = reshape(task_cost.', 1, []);
                    
                    obj_driver = 0;
                    LP_time_driver = zeros(num_node^2, 1);
                    
                    filename = strcat("result/result_proposed/num_driver_", num2str(sum_driver), "/node_", num2str(num_node), "/depot_", num2str(num_depot), "/net_", num2str(net_num), "/allocation/allocation_", num2str(exp_num));
                    
                    driver_allocation_filename = strcat(filename, "/allocation_od");
                    task_allocation_filename = strcat(filename, "/allocation_rs");
                    mkdir(driver_allocation_filename);
                    mkdir(task_allocation_filename);
                    
                    driver_allocation_filename = strcat(filename, "/allocation_od/driver_allocation_", num2str(sum_driver), "_", num2str(num_node), "_", num2str(num_depot), "_", num2str(net_num), "_", num2str(exp_num));
                    task_allocation_filename = strcat(filename, "/allocation_rs/task_allocation_", num2str(sum_driver), "_", num2str(num_node), "_", num2str(num_depot), "_", num2str(net_num), "_", num2str(exp_num));
                    
                    
                    clear filename;
                    
                    for od_index = 1:num_node^2
                        
                        driver_index_first = sum(num_driver(1:od_index-1)) + 1;
                        driver_index_last = sum(num_driver(1:od_index));
                        
                        
                        %odドライバ�?�に配�?されたタスク数を取�?
                        task_od = driver_allocation(od_index, 1:num_depot*num_node+1);
                        
                        %配�?されたタスクが非0である要�?を抜き�?��?
                        [row_task, col_task, v_task] = find(task_od);
                        
                        %ドライバ�?�の知覚費用を取�?
                        driver_cost_od = driver_cost(driver_index_first:driver_index_last, col_task);
                        
                        
                        %目�?関数の係数行�??
                        %-------------------------------------------------------------------------------------------------------------------------------------------
                        f = reshape(driver_cost_od.', 1, []);
                        
                        %等式制�?条件の右辺
                        %---------------------------------------------------------------------------------------------------------------------------------------------------
                        b_1 = v_task;
                        b_2 = ones(1, num_driver(od_index));
                        beq = [b_1, b_2];
                        
                        clear b_1;
                        clear b_2;
                        
                        %等式制�?条件の係数行�??
                        %--------------------------------------------------------------------------------------------------------------------------------------------------------
                        data_1 = ones(1, num_driver(od_index)*numel(col_task));
                        data_2 = ones(1, num_driver(od_index)*numel(col_task));
                        data_eq = [data_1, data_2];
                        
                        clear data_1;
                        clear data_2;
                        
                        row_1 = zeros(1, num_driver(od_index)*numel(col_task));
                        for task_index = 1:numel(col_task)
                            for j = 1:num_driver(od_index)
                                row_1((task_index-1)*num_driver(od_index) + j) = task_index;
                            end
                        end
                        row_2 = zeros(1, num_driver(od_index)*numel(col_task));
                        for i = 1:num_driver(od_index)
                            for task_index = 1:numel(col_task)
                                row_2((i-1)*numel(col_task) + task_index) = numel(col_task) + i;
                            end
                        end
                        row_eq = [row_1, row_2];
                        
                        clear row_1;
                        clear row_2;
                        
                        col_1 = zeros(1, num_driver(od_index)*numel(col_task));
                        for task_index = 1:numel(col_task)
                            for i = 1:num_driver(od_index)
                                col_1((task_index-1)*num_driver(od_index) + i) = (i-1)*(numel(col_task)) + task_index;
                            end
                        end
                        col_2 = 1:num_driver(od_index)*numel(col_task);
                        col_eq = [col_1, col_2];
                        
                        clear col_1;
                        clear col_2;
                        
                        Aeq = sparse(row_eq, col_eq, data_eq);
                        clear row_eq;
                        clear col_eq;
                        clear data_eq;
                        
                        LP_time = cputime;
                        
                        [x, fval] = linprog(f, [], [], Aeq, beq, zeros(1, num_driver(od_index)*(numel(col_task))), ones(1, num_driver(od_index)*numel(col_task)));
                         
                        LP_time_driver(od_index) = cputime - LP_time;
                        obj_driver = obj_driver + fval;
                        
                        %結果の書き�?��?
                        row_od = zeros(1, numel(x));
                        for i = 1:num_driver(od_index)
                            for j = 1:numel(col_task)
                                row_od(numel(col_task)*(i-1) + j) = i;
                            end
                        end
                        col_od = zeros(1, numel(x));
                        for i = 1:num_driver(od_index)
                            for j = 1:numel(col_task)
                                col_od(numel(col_task)*(i-1) + j) = col_task(j);
                            end
                        end
                        atomic_allocation_od = sparse(row_od, col_od, x, num_driver(od_index), num_depot*num_node+1);
                        
                        clear row_od;
                        clear col_od;
                        
                        driver_allocation_filename_od = strcat(driver_allocation_filename, "_od_index_", num2str(od_index-1), ".csv");
                        
                        writematrix(atomic_allocation_od, driver_allocation_filename_od);
                        
                        clear atomic_allocation_od;
                        clear driver_allocation_filename_od;
                       
                        
                        clear x;
                        clear fval;
                        
                        clear f;
                        clear Aeq;
                        clear beq;
                        
                        clear row_task;
                        clear col_task;
                        clear v_task;
                        
                        clear driver_cost_od
                        
                        clear driver_index_first;
                        clear driver_index_last;
                        
                    end
                    
                    %LP_time_driver
                    
                    "Complete Driver"
                    
                    clear driver_od;
                    
                    LP_time_task = zeros(num_depot*num_node, 1);
                    obj_task = 0;
                    
                    for rs_index = 1:num_depot*num_node
                        
                        task_index_first = sum(num_task(1:rs_index-1)) + 1;
                        task_index_last = sum(num_task(1:rs_index));
                        
                        %rsタスク�?有�??に配�?された実行されるタスク数(実行されな�?タスク数)を取�?
                        task_rs = [task_allocation(rs_index, 2), task_allocation(rs_index, 1)];
                        
                        %タスク�?有�??の知覚費用を取�?
                        task_cost_rs = task_cost(task_index_first:task_index_last, 1:2);
                        
                        %目�?関数の係数行�??
                        %-------------------------------------------------------------------------------------------------------------------------------------------
                        f = reshape(task_cost_rs.', 1, []);
                        
                        %等式制�?条件の右辺
                        %---------------------------------------------------------------------------------------------------------------------------------------------------
                        b_1 = task_rs;
                        b_2 = ones(1, num_task(rs_index));
                        beq = [b_1, b_2];
                        
                        clear b_1;
                        clear b_2;
                        
                        %等式制�?条件の係数行�??
                        %--------------------------------------------------------------------------------------------------------------------------------------------------------
                        data_1 = ones(1, num_task(rs_index)*2);
                        data_2 = ones(1, num_task(rs_index)*2);
                        data_eq = [data_1, data_2];
                        
                        clear data_1;
                        clear data_2;
                        
                        row_1 = zeros(1, num_task(rs_index)*2);
                        for index = 1:2
                            for j = 1:num_task(rs_index)
                                row_1((index-1)*num_task(rs_index) + j) = index;
                            end
                        end
                        row_2 = zeros(1, num_task(rs_index)*2);
                        for i = 1:num_task(rs_index)
                            for index = 1:2
                                row_2((i-1)*2 + index) = 2 + i;
                            end
                        end
                        row_eq = [row_1, row_2];
                        
                        clear row_1;
                        clear row_2;
                        
                        col_1 = zeros(1, num_task(rs_index)*2);
                        for index = 1:2
                            for i = 1:num_task(rs_index)
                                col_1((index-1)*num_task(rs_index) + i) = (i-1)*2 + index;
                            end
                        end
                        col_2 = 1:num_task(rs_index)*2;
                        col_eq = [col_1, col_2];
                        
                        clear col_1;
                        clear col_2;
                        
                        Aeq = sparse(row_eq, col_eq, data_eq);
                        clear row_eq;
                        clear col_eq;
                        clear data_eq;
                        
                        LP_time = cputime;
                        
                        [x, fval] = linprog(f, [], [], Aeq, beq, zeros(1, num_task(rs_index)*2), ones(1, num_task(rs_index)*2));
                           
                        LP_time_task(rs_index) = cputime - LP_time;
                        obj_task = obj_task + fval;
                        
                        %結果の書き�?��?
                        atomic_allocation_rs = reshape(x, 2, []).';
                        task_allocation_filename_rs = strcat(task_allocation_filename, "_rs_index_", num2str(rs_index-1), ".csv");
                        
                        writematrix(atomic_allocation_rs, task_allocation_filename_rs);
                        
                        clear atomic_allocation_rs;
                        clear task_allocation_filename_rs;
                        
                        clear x;
                        clear fval;
                        
                        clear f;
                        clear Aeq;
                        clear beq;
                        
                        clear task_cost_od
                        
                        clear task_index_first;
                        clear task_index_last;
                        
                    end
                    
                    clear driver_allocation_filename;
                    clear task_allocation_filename;
                    
                    %LP_time_task
                    
                    "Complete Task"
                    
                    clear task_rs;
                    
                    obj_driver;
                    obj_task;
                    obj = obj_driver + obj_task;
                    
                    proposal_obj_array(net_num+1, exp_num+1) = obj;
                    
                    filename = strcat("result/result_proposed/num_driver_", num2str(sum_driver), "/node_", num2str(num_node), "/depot_", num2str(num_depot), "/time/net_", num2str(net_num), "/exp_", num2str(exp_num));
                    mkdir(filename);
                    
                    driver_time_filename = strcat(filename, "/LP_time_driver.csv");
                    task_time_filename = strcat(filename, "/LP_time_task.csv");
                    
                    writematrix(LP_time_driver, driver_time_filename);
                    writematrix(LP_time_task, task_time_filename);
                    
                    clear driver_time_filename;
                    clear task_time_filename;
                    
                    LP_time = [LP_time_driver; LP_time_task];
                    clear LP_time_driver;
                    clear LP_time_task;
                    
                    max_LP_time_array(net_num+1, exp_num+1) = max(LP_time);
                    sum_LP_time_array(net_num+1, exp_num+1) = sum(LP_time);
                    
                    clear LP_time;
                    
                    clear obj_driver;
                    clear obj_task;
                    clear obj;
                    
            
                    writematrix(max_LP_time_array, max_time_filename);
                    writematrix(sum_LP_time_array, sum_time_filename);
                    writematrix(proposal_obj_array, obj_filename);
                    
                end
                
            end
            
            clear foldername;
            
            
            
            clear max_LP_time_array;
            clear sum_LP_time_array;
            clear proposal_obj_array;
            clear filename;
            clear max_time_filename;
            clear sum_time_filename;
            clear obj_filename;
            
            
        end
    end 
end

clear num_driver_array;
clear num_node_array;
clear num_depot_array;

clear num_networks;
clear num_experience;