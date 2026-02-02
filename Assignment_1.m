function [R, rnk, E_list, E_total] = my_rref(A)
    [m, n] = size(A);
    R = A;
    E_list = {};
    E_total = eye(m);
    cnt = 1;

    pr = 1;
    pcols = [];
    for c = 1:n
        if pr > m
            break;
        end

        [~, mi] = max(abs(R(pr:m, c)));
        mi = mi + pr - 1;

        if abs(R(mi, c)) < 1e-10
            continue;
        end
        if mi ~= pr
            E = eye(m);
            E([pr, mi], :) = E([mi, pr], :);
            E_list{cnt} = E;
            E_total = E * E_total;
            R = E * R;
            cnt = cnt + 1;
        end
        if abs(R(pr, c) - 1) > 1e-10
            sf = R(pr, c);
            E = eye(m);
            E(pr, pr) = 1 / sf;
            E_list{cnt} = E;
            E_total = E * E_total;
            R = E * R;
            cnt = cnt + 1;
        end
        for r = pr+1:m
            if abs(R(r, c)) > 1e-10
                f = R(r, c);
                E = eye(m);
                E(r, pr) = -f;
                E_list{cnt} = E;
                E_total = E * E_total;
                R = E * R;
                cnt = cnt + 1;
            end
        end

        pcols = [pcols, c];
        pr = pr + 1;
    end
    for i = length(pcols):-1:1
        pc = pcols(i);
        pr = i;

        for r = 1:pr-1
            if abs(R(r, pc)) > 1e-10
                f = R(r, pc);
                E = eye(m);
                E(r, pr) = -f;
                E_list{cnt} = E;
                E_total = E * E_total;
                R = E * R;
                cnt = cnt + 1;
            end
        end
    end

    R(abs(R) < 1e-10) = 0;
    E_total(abs(E_total) < 1e-10) = 0;
    rnk = sum(any(abs(R) > 1e-10, 2));
    fprintf('\nRREF:\n');
    disp(R);

    fprintf('Rank = %d\n', rnk);

    fprintf('\nFinal transformation matrix E_total:\n');
    disp(E_total);
end

A = [2, 1, -1, 8;
     -3, -1, 2, -11;
     -2, 1, 2, -3];
fprintf('modify A in the code if needed');
tic;
[R, rnk, E_list, E_total] = my_rref(A);
t1 = toc;

tic;
R2 = rref(A);
t2 = toc;

fprintf('My time: %.6f sec\n', t1);
fprintf('Built-in time: %.6f sec\n', t2);
disp('Built-in result:');
disp(R2);


function draw_system(A, b, name)
    fprintf('\nSystem %s\n', name);
    [m, n] = size(A);
    Aug = [A, b];
    [R_aug, ~, ~] = my_rref(Aug);    
    disp('RREF of augmented:');
    disp(R_aug);    
    rA = sum(any(abs(R_aug(:,1:n)) > 1e-10, 2));
    rAb = sum(any(abs(R_aug) > 1e-10, 2));
    fprintf('rank(A) = %d, rank([A|b]) = %d\n', rA, rAb);
    
    if rA < rAb
        sol = 'No solution';
    elseif rA == n
        sol = 'Unique solution';
    else
        sol = sprintf('%d free variables', n - rA);
    end    
    fprintf('Solution: %s\n', sol);
    if n == 3
        figure('Name', name);
        subplot(1, 2, 1);
        hold on; grid on;
        title(['Planes - ' name]);
        xlabel('x'); ylabel('y'); zlabel('z');
        view(45, 30);        
        [xg, yg] = meshgrid(-5:0.5:5, -5:0.5:5);
        cols = ['r', 'g', 'b', 'm'];        
        for i = 1:m
            if abs(A(i, 3)) > 1e-10
                zg = (b(i) - A(i,1)*xg - A(i,2)*yg) / A(i,3);
            elseif abs(A(i, 2)) > 1e-10
                yg2 = (b(i) - A(i,1)*xg) / A(i,2);
                zg = zeros(size(xg));
                surf(xg, yg2, zg, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', cols(i));
                continue;
            else
                continue;
            end
            surf(xg, yg, zg, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', cols(i));
        end
        
        subplot(1, 2, 2);
        hold on; grid on;
        title(['Vectors - ' name]);
        xlabel('x'); ylabel('y'); zlabel('z');
        view(45, 30);
        
        for i = 1:n
            quiver3(0, 0, 0, A(1,i), A(2,i), A(3,i), 'LineWidth', 2, 'Color', cols(i));
        end
        
        quiver3(0, 0, 0, b(1), b(2), b(3), 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');
        
        legend('a1', 'a2', 'a3', 'b');
    end
end
Aa = [1, 1, 1; 1, 2, 3; 0, 1, 2];
ba = [2; 1; 0];
draw_system(Aa, ba, 'a');

Ab = [1, 1, 1; 1, 2, 1; 2, 3, 2];
bb = [2; 3; 5];
draw_system(Ab, bb, 'b');

Ac = [1, 1, 1; 1, 2, 1; 2, 3, 2];
bc = [2; 3; 9];
draw_system(Ac, bc, 'c');

Ad = [1, 1, 1; 1, 2, 1; 2, 9, 5];
bd = [2; 3; 12];
draw_system(Ad, bd, 'd');


function x = backsub(U, c)
    n = length(c);
    x = zeros(n, 1);
    
    for i = n:-1:1
        pc = find(abs(U(i, :)) > 1e-10, 1);
        if isempty(pc)
            continue;
        end
        
        x(pc) = (c(i) - U(i, pc+1:n) * x(pc+1:n)) / U(i, pc);
    end
end

function x = solve_system(A, b)
    [~, n] = size(A);
    Aug = [A, b];
    
    [R_aug, ~, ~] = my_rref(Aug);
    
    U = R_aug(:, 1:n);
    c = R_aug(:, n+1);
    
    rA = sum(any(abs(U) > 1e-10, 2));
    rAug = sum(any(abs(R_aug) > 1e-10, 2));
    
    if rA < rAug
        error('inconsistent');
    end
    
    x = backsub(U, c);
end

% test
At = [1, 1, 1; 1, 2, 1; 2, 9, 5];
bt = [2; 3; 12];

disp('Test A:');
disp(At);
disp('Test b:');
disp(bt);

tic;
x1 = solve_system(At, bt);
t1 = toc;
disp('My solution:');
disp(x1);
fprintf('Time: %.6f sec\n', t1);

tic;
x2 = At\bt;
t2 = toc;
disp('inv(A)*b:');
disp(x2);
fprintf('Time: %.6f sec\n', t2);

tic;
x3 = At \ bt;
t3 = toc;
disp('A\b:');
disp(x3);
fprintf('Time: %.6f sec\n', t3);

disp('Residuals:');
disp(At * x1 - bt);
disp(At * x2 - bt);
disp(At * x3 - bt);

fprintf('My time: %.6f\n', t1);
fprintf('inv time: %.6f (%.2fx)\n', t2, t1/t2);
fprintf('backslash time: %.6f (%.2fx)\n', t3, t1/t3);