
% 個人のコスト行列(num_driver*num_rs)の読み込み
indivisual_cost = readmatrix("input/indivisual_cost.csv");
num_shippers = readmatrix("input/num_shippers.csv");
num_drivers = readmatrix("input/num_drivers.csv");


Aeq = make_Aeq(indivisual_cost);
    clear Aeq_row;
    clear Aeq_col;
    clear Aeq_value;
beq = ones(1, sz(1));

A = make_A(indivisual_cost);
    clear A_row;
    clear A_col;
    clear A_value;
b = reshape((-1)*num_shippers.', 1, []);

lb = zeros(1, sz(1)*sz(2));
ub = ones(1, sz(1)*sz(2));

f = reshape(indivisual_cost.', 1, []);

[x, fval] = linprog(f, A, b, Aeq, beq, lb, ub)
x = reshape(x, sz(2), sz(1)).';

mkdir("output");
writematrix(x, "output/atomic_solution.csv");
writematrix(fval, "output/func_value.csv");




function Aeq = make_Aeq(indivisual_cost)
    sz = size(indivisual_cost);
    row = ones(1, sz(2));
    Aeq_row = row;
    for driver = 2:sz(1)
        add_row = row*driver;
        Aeq_row = cat(2, Aeq_row, add_row);
    end
    clear row;
    clear add_row;
    Aeq_col = 1:sz(1)*sz(2);
    Aeq_value = ones(1, sz(1)*sz(2));
    Aeq = sparse(Aeq_row, Aeq_col, Aeq_value); 
end

function A = make_A(indivisual_cost)
    sz = size(indivisual_cost);
    row = 1:sz(2);
    A_row = row;
    for driver = 2:sz(1)
        A_row = cat(2, A_row, row);
    end
    clear row;
    A_col = 1:sz(1)*sz(2);
    A_value = ones(1, sz(1)*sz(2))*(-1);
    A = sparse(A_row, A_col, A_value);
end