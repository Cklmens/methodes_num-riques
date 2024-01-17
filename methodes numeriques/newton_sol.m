function x = newton_sol(B, b, g, x0, tolerance, max_iterations)
    % Initialisation
    x = x0;
    iteration = 0;
    
    while true
        % Calcul de F(x)
        F = min(B*x - b, x - g);
      
        % Vérification de la condition d'arrêt
        if norm(F, inf) < tolerance || iteration >= max_iterations
            break;
        end
        
        % Calcul de la dérivée F'_0(x)
        F_prime = zeros(size(B));
        for i = 1:size(B, 1)
            tmp=B*x - b;
            temp=(x - g);
            if tmp(i) <= temp(i)
                F_prime(i, :) = B(i, :);
            else
                F_prime(i, i) = 1;
            end
        end
        
        % Mise à jour de x selon la méthode de Newton
        x = x - inv(F_prime)*F ;
        
        % Incrémentation du nombre d'itérations
        iteration = iteration + 1;
    end
end