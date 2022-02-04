per_driver_array = [128, 256, 512];
num_node_array = [25];
num_depot_array = [8];

num_networks = 2;
num_experience = 2;

    
for num_node = num_node_array
        
    for num_depot = num_depot_array
               
        for per_driver = per_driver_array
    
            sum_driver = per_driver * num_node * num_node;
            sum_task = sum_driver;
            
            filename = strcat("../../data/nn_", num2str(num_node), "/nv_", num2str(per_driver), "/nd_", num2str(num_depot));
            time_filename = strcat(filename, "/time_summary/original/original_time.csv");
            obj_filename = strcat(filename, "/obj_summary/original/original_obj.csv");
            
            foldername = strcat(filename, "/time_summary");
            mkdir(foldername);
            foldername = strcat(foldername, "/original");
            mkdir(foldername);
            
            foldername = strcat(filename, "/obj_summary");
            mkdir(foldername);
            foldername = strcat(foldername, "/original");
            mkdir(foldername);
            
            clear foldername;
            
            LP_time_array = zeros(num_networks, num_experience);
            LP_obj_array = zeros(num_networks, num_experience);
            
            for net_num = 0:num_networks-1
                for exp_num = 0:num_experience-1
   
                    filename = strcat("../../data/nn_", num2str(num_node), "/nv_", num2str(per_driver), "/nd_", num2str(num_depot), "/net_", num2str(net_num), "/exp_", num2str(exp_num));
                    driver_cost_atomic_filename = strcat(filename, "/input/cost_atomic/driver_cost.csv");
                    task_cost_atomic_filename = strcat(filename, "/input/cost_atomic/shipper_cost.csv");
                    num_driver_filename = strcat(filename, "/input/num_drivers.csv");
                    num_task_filename = strcat(filename, "/input/num_shippers.csv");
                    
                    driver_cost = readmatrix(driver_cost_atomic_filename);
                    task_cost = readmatrix(task_cost_atomic_filename);
                    num_driver = round(readmatrix(num_driver_filename));
                    num_task = round(readmatrix(num_task_filename));
                    
                    
                    num_driver = reshape(num_driver.', 1, []);
                    num_task = reshape(num_task.', 1, []);
                    
                    clear filename;
                    clear driver_cost_atomic_filename;
                    clear task_cost_atomic_filename;
                    clear demand_filename;
                    
                    driver_cost = reshape(driver_cost.', 1, []);
                    task_cost = reshape(task_cost.', 1, []);
                    
                    %目�?関数の係数行�??
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                    f = [driver_cost, task_cost];
                    
                    %不等式制�?条件の右辺
                    %------------------------------------------------------------------------------
                    b = zeros(1, num_depot*num_node);
                    
                    %等式制�?条件の右辺
                    %-------------------------------------------------------------------------------
                    beq = ones(1, sum_driver + sum_task);
                    
                    
                    %不等式制�?条件の係数行�??
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    data_1 = ones(1, sum_driver*num_depot*num_node);
                    data_1 = data_1.*(-1);
                    data_2 = ones(1, sum_task);
                    data = [data_1, data_2];
                        
                    clear data_1;
                    clear data_2;
                    
                    row_1 = zeros(1, sum_driver*num_depot*num_node);
                    for i = 1:num_depot*num_node
                        for j = 1:sum_driver
                            row_1((i-1)*sum_driver + j) = i;
                        end
                    end
                    row_2 = zeros(1, sum_task);
                    tasks = 0;
                    for rs_index = 1:num_depot*num_node
                        for j = 1:num_task(rs_index)
                            row_2(tasks + j) = rs_index;
                        end
                        tasks = tasks + num_task(rs_index);
                    end
                    row = [row_1, row_2];
                    
                    clear row_1;
                    clear row_2;
                    clear tasks;
                   
                    col_1 = zeros(1, sum_driver*num_depot*num_node);
                    for rs_index = 1:num_depot*num_node
                        for i = 1:sum_driver
                            col_1((rs_index-1)*sum_driver + i) = (i-1)*(num_depot*num_node+1) + rs_index;
                        end
                    end
                    col_2 = sum_driver*(num_depot*num_node+1);
                    for i = 1:sum_task
                        col_2(i) = sum_driver*(num_depot*num_node+1) + 2*i;
                    end
                    col = [col_1, col_2];
                    
                    clear col_1;
                    clear col_2;
                    
                    A = sparse(row, col, data);
                    
                    
                    
                    %等式制�?条件の係数行�??
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    data_eq = ones(1, sum_driver*(num_depot*num_node + 1) + sum_task*2);
                    
                    
                    row_eq = zeros(1, sum_driver*(num_depot*num_node + 1) + sum_task*2);
                    for i = 1:sum_driver
                        for rs_index = 1:num_depot*num_node+1
                            row_eq((i-1)*(num_depot*num_node+1) + rs_index) = i;
                        end
                    end
                    for j = 1:sum_task
                        row_eq(sum_driver*(num_depot*num_node+1) + 2*(j-1) + 1) = sum_driver + j;
                        row_eq(sum_driver*(num_depot*num_node+1) + 2*(j-1) + 2) = sum_driver + j;
                    end                     
                    
                    col_eq = 1:sum_driver*(num_depot*num_node + 1) + sum_task*2;
                    
                    %不等号ので計算を実行する�?�合には以下�?�まとまりをコメントアウトすること
                    %---------------------------------------------------------------------------------------------------------------------------------------------------
                    
                    row_plus = ones(1, sum_driver*(num_depot*num_node + 1) + sum_task*2).*(num_depot*num_node);
                    
                    row_eq = row_eq + row_plus;
                    
                    clear row_plus;
                    
                    row_eq = [row, row_eq];
                    col_eq = [col, col_eq];
                    data_eq = [data, data_eq];
                    beq = [b, beq];
                    
                    %----------------------------------------------------------------------------------------------------------------------------------------------------
                    
                    Aeq = sparse(row_eq, col_eq, data_eq);
                    %size(Aeq)
                    %size(beq)
                    
                    clear data_eq;
                    clear row_eq;
                    clear col_eq;
                    
                    clear data;
                    clear row;
                    clear col;
                    
                    
                    
                    
                    
                    LP_time = cputime;
                    
                    %[x, fval] = linprog(f, A, b, Aeq, beq, zeros(1, sum_driver*(num_depot*num_node+1) + sum_task*2), ones(1, sum_driver*(num_depot*num_node+1) + sum_task*2));
                    [x, fval] = linprog(f, [], [], Aeq, beq, zeros(1, sum_driver*(num_depot*num_node+1) + sum_task*2), ones(1, sum_driver*(num_depot*num_node+1) + sum_task*2));
                
                    LP_time_array(net_num+1, exp_num+1) = cputime - LP_time;
                    
                    
                    
                    
                    clear f;
                    clear A;
                    clear b;
                    clear Aeq;
                    clear beq;
                    clear cputime;
                    clear LP_time;
                    
                    
                    atomic_allocation_driver = x(1:sum_driver*(num_depot*num_node+1));
                    atomic_allocation_task = x(sum_driver*(num_depot*num_node+1)+1:sum_driver*(num_depot*num_node+1)+sum_task*2);
                    
                    atomic_allocation_driver = reshape(atomic_allocation_driver, num_depot*num_node+1, []).';
                    atomic_allocation_task = reshape(atomic_allocation_task, 2, []).';
                    
                    foldername = strcat("../../data/nn_", num2str(num_node), "/nv_", num2str(per_driver), "/nd_", num2str(num_depot), "/net_", num2str(net_num), "/exp_", num2str(exp_num), "/output");
                    mkdir(foldername);
                    foldername = strcat(foldername, "/original");
                    mkdir(foldername);
                    foldername = strcat(foldername, "/sol");
                    mkdir(foldername);
                    foldername = strcat(foldername, "/sol_allocation");
                    mkdir(foldername);
                    
                    driver_allocation_filename = strcat(foldername, "/driver_allocation.csv");
                    task_allocation_filename = strcat(foldername, "/shipper_allocation.csv");
                    
                    
                    LP_obj_array(net_num+1, exp_num+1) = fval;
                    clear fval;
                    
                    writematrix(atomic_allocation_driver, driver_allocation_filename);
                    writematrix(atomic_allocation_task, task_allocation_filename);
                    
                    clear atomic_allocation_driver;
                    clear atomic_allocation_task;
                    clear driver_allocation_filename;
                    clear task_allocation_filename;
                    
                    
                    writematrix(LP_time_array, time_filename);
                    writematrix(LP_obj_array, obj_filename);
                    
                    
                    
                end
            end
            
            
            clear foldername;
            
            
            
            
            clear LP_time_array;
            clear LP_obj_array;
            clear filename;
            clear time_filename;
            clear obj_filename;
            
        end
    end
end