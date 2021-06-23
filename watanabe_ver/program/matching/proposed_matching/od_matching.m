% 個人のコスト行列(num_driver*num_rs)の読み込み
indivisual_cost = readmatrix("input/indivisual_cost.csv");
num_drivers = readmatrix("input/num_drivers.csv");
task_allocation = readmatrix("input/int_allocation.csv");


sz_dv = size(num_drivers);
sz = size(indivisual_cost);
first_driver_index = 1;

time = zeros(sz_dv(1), sz_dv(2));

sum_obj = 0;

mkdir("output");

whole_task_allocation = zeros(1, sz(2));


for o = 1:sz_dv(1)
    for d = 1:sz_dv(2)
        
        if num_drivers(o, d) == 0
            continue;
        end
        
        od = sz_dv(2)*(o-1) + d;
        n = task_allocation(od, :);
        
        beq = make_beq(num_drivers(o, d), n);
        
        Aeq = make_Aeq(num_drivers(o, d), sz(2));
            clear Aeq_row;
            clear Aeq_col;
            clear Aeq_value;
        
        lb = zeros(1, num_drivers(o, d)*sz(2));
        ub = ones(1, num_drivers(o, d)*sz(2));
        
        atomic_cost_od = indivisual_cost(first_driver_index:first_driver_index+num_drivers(o, d)-1, :);
        f = reshape(atomic_cost_od.', 1, []);
        
        size(f);
        size(Aeq);
        
        
        tic

        [x, fval] = linprog(f, [], [], Aeq, beq, lb, ub);
        x = reshape(x, sz(2), num_drivers(o, d)).';
        
        whole_task_allocation = cat(1, whole_task_allocation, x);
                
        time(o, d) = toc;
        
        
        
        mkdir("output/od(" + od + ")");
        writematrix(x, "output/od(" + od + ")/atomic_solution.csv");
        writematrix(fval, "output/od(" + od + ")/func_value.csv");
        
        sum_obj = sum_obj + fval;
        
        clear Aeq;
        clear beq;
        clear lb;
        clear ub;
        clear f;
        
        first_driver_index = first_driver_index + num_drivers(o, d);
        
        str = "complete od = " + o + "×" + d
        
            
    end
end

whole_task_allocation(1, :) = [];

writematrix(whole_task_allocation, "output/whole_task_allocation.csv");
writematrix(sum_obj, "output/sum_func_value.csv");
writematrix(time, "output/time.csv");



function beq = make_beq(num_driver, task_allocation)
    
    beq = ones(1, num_driver);
    beq = cat(2, beq, task_allocation);

end

function Aeq = make_Aeq(num_driver, num_rs)

    row = ones(1, num_rs);
    Aeq_row = row;
    for driver = 2:num_driver
        add_row = row*driver;
        Aeq_row = cat(2, Aeq_row, add_row);
    end
    
    row = num_driver + 1:num_driver + num_rs;
    for driver = 1:num_driver
        Aeq_row = cat(2, Aeq_row, row);
    end
    
    size(Aeq_row);
    
    clear row;

    Aeq_col = 1:num_driver*num_rs;
    Aeq_col = cat(2, Aeq_col, Aeq_col);
    
    size(Aeq_col);
    
    Aeq_value = ones(1, num_driver*num_rs);
    Aeq_value = cat(2, Aeq_value,  ones(1, num_driver*num_rs));
    
    Aeq = sparse(Aeq_row, Aeq_col, Aeq_value);

end