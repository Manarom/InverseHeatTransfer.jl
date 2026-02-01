classdef AbstractPolynomial<handle
    %
    properties (Constant)
        left = -1
        right = 1
    end 
    properties (SetAccess = protected)
        coeffs % coefficients of current polynomial
        xmin % min of x initial
        xmax % max of x_initial
        x_normalized % normalized to -1...1 coordinates used to create the vandermatrix
        x_initial = [] % initial coordinates before the normalization, used to create the vandermatrix 
    end
    
    properties (SetAccess = protected)
        vandermatrix % Vandermonde matrix
        R % 
        Q % QR factorization of the vandermonde matrix used speedup the fitting
    end
    properties (Dependent)
        degree % polynomial degree
        N % polynomial coeffs number
        npoints % number of coordinates
    end
    methods(Abstract)
        monomial(obj,i,x) % evaluates monomial
        monomial_derivative(obj,i,x) % derivative of the monomial
        derivative_conversion(obj,coeffs) % function to convert the coefficients of the initial polynomial to the coefficients 
        % of polynomial deriavative, coefficients of the derivative polynomial are in the same bassis of degree - 1
    end
    methods
        function val = get.degree(obj)
            val = obj.N - 1;
        end
        function val = get.N(obj)
            val = numel(obj.coeffs);
        end
        function val= get.npoints(obj)
            val = numel(obj.x_normalized);
        end
        function set_coeffs(obj,co)
            % sets new coefficients 
            n = numel(co);
            if n==0
                return
            end
            if n ~= obj.N
                fill_vander_matrix(obj,obj.x_initial,n - 1)
            end
            obj.coeffs = co(:);
        end
        function obj = AbstractPolynomial(x,degree)
            obj.coeffs = ones(degree + 1,1);
            fill_vander_matrix(obj,x,degree)
        end
        function cpObj = deep_copy(ply,new_coordinates,new_coeffs)
            arguments
                ply
                new_coordinates =[] % new coordinates vector to evaluate the vandermatrix, if empty copies 
                % the one from ply
                new_coeffs = []
            end
            cls = class(ply);
            if ~isempty(new_coordinates)
                x = new_coordinates;
            else
                x = ply.x_initial;
            end
            if ~isempty(new_coeffs)
                cs = ply.coeffs;
            else
                cs = new_coeffs;
            end
            cpObj = feval(cls, x, length(cs)-1);  % Calls subclass constructor with no args
            set_coeffs(cpObj, cs);
        end
        function  val = derivative_degree(obj)
            val = obj.degree - 1;
        end
        function stat = distance(p1,p2,x,y)
            % evaluates distance between two polynomials 
            % p1 - first polynomial
            % p2 - second polynomial
            % x - independent variable if it is empty methods takes it from
            % the initial polynomial
            arguments
                p1
                p2 AbstractPolynomial
                x = []
                y = []
            end
            is_y_empty = isempty(y);
            is_x_empty = isempty(x);
            if ~is_x_empty && ~is_y_empty
                assert(length(x)==length(y), "x and y must be of the same size")
            end
            d1 = p1.degree;
            d2 = p2.degree;
            if is_x_empty
                x = p1.x_initial; 
            else
                fill_vander_matrix(p1,x,d1) % filling first poly
            end
            c1 = p1.coeffs;
            y1 = poly_eval_unnorm(p1,c1,x);
            if is_y_empty
                y1 = y;
            end
            r1 = y - y1;
            fill_vander_matrix(p2,x,d2) % filling second poly
            c2 = poly_fit(p2,y);
            y2 = poly_eval_unnorm(p2,c2,x);
            r2 = y - y2; 
            r12 = y1 - y2;
            stat = struct("r12",r12,"r1",r1,"r2",r2,"c1",c1,"c2",c2);
        end
        function y = poly_eval_unnorm(obj,coeffs,x)
            % evaluates polynomial for unnormalized x!
            % if x lies outside the polynomial basis range 
            % returns error 
            x = x(:);
            M = numel(x);
            n_init = obj.N;
            n_cur = length(coeffs);
            
            if (n_init == n_cur) && isequal(obj.x_initial,x)  % if the degree of polynomial does not change
                 y = obj.vandermatrix*coeffs(:);
                return
            end
            [a,b] = AbstractPolynomial.scalers();
            [xm,xx] = obj.x_scalers();
            n = n_cur - 1;
            y = zeros(M,1);
            x = a + (b - a)*(x - xm)/(xx - xm);
            flag = (x <= b) & (x >= a); 
            x(x < a) = a;
            x(x > b) = b;
            if all(flag)
                for i = 0:n
                    y = y + monomial(obj,n,i,x)*coeffs(i + 1); % just evaluating monomials without changing the obj itself
                end
                return
            end
            for i = 0:n
                y(flag) = y(flag) + monomial(obj,n,i,x(flag))*coeffs(i + 1); % just evaluating monomials without changing the obj itself
            end
            y(~flag) = interp1(x(flag),y(flag),x(~flag),"pchip","extrap");
        end
        function poly = poly_der(obj)
            % returns the polynomial of the same type 
            % which is the derivative of the initial poynomial with respect to the independent 
            % variable
            ActualType = class(obj); % getting actual type of the object 
            poly = feval(ActualType, obj.x_initial, obj.degree - 1); 
            poly.coeffs = derivative_conversion(poly,obj.coeffs);

        end
        function J = jacobian(obj,x)
            % evaluates jacobian for x within the range of polynomial basis
            % otherwise returns an error
            if isequal(x,obj.x_initial)
                J = obj.vandermatrix;
            else
                if all(x<=obj.xmax) && all(x>=obj.xmin)
                    normalized_x  = AbstractPolynomial.normalize_x(x(:));
                    eval_vandermatrix(obj,normalized_x, obj.degree);
                else
                    error("to evaluate the jacobian x values must be within the xmin xmax range for polynomial")
                end
            end
        end
        function [mn,mx] = x_scalers(obj)
            mn = obj.xmin;
            mx = obj.xmax;
        end
        function [coefs,V] = poly_fit(obj,y,options)
            arguments
                obj
                y
                options.x = []
            end
            if isempty(options.x) || isequal(options.x,obj.x_initial)
                coefs = obj.R\(transpose(obj.Q)*y);
                obj.coeffs = coefs;
                if nargout > 1
                    V = obj.vandermatrix;
                end
                return
            end
            y = y(:);
            x = normalize_x(options.x(:));
            V = eval_vandermatrix(obj, x, obj.degree);
            coefs = V\y;
        end
        function fun = as_fun(obj)
            % wraps obj as a function, which modifies both x and P 
            fun = @(x,coeffs) poly_eval_unnorm(obj,coeffs,x);     % Cp_parametric_fun = @(T,P)polyval(P,T);
        end
        function ax = plot(obj,options)
            arguments 
                obj
                options.ax = []
                options.basis (1,1) logical = false % if true  - plot basis function
                options.value (1,1) logical = true % if true plot function value
                options.hold {mustBeMember(options.hold, ["on", "off"])} = "off"
            end
            
            if isempty(options.ax)
                ax = get_next_ax();
            else
                ax = options.ax;
            end
            if ~options.basis && ~options.value
                return
            end
            is_hold = options.hold == "on";
            if is_hold
                hold(ax,"on")
            end
            if options.basis
                plot(ax,obj.x_initial,obj.vandermatrix);
            end
            if options.value
                if options.basis && ~is_hold
                    hold(ax,"on")
                end
                plot(ax,obj.x_initial,obj.vandermatrix*obj.coeffs);
            end
            hold(ax,"off")
        end
    end
    methods (Access = protected)
        function fill_x(obj,x)
            % filling independent variable from
            obj.x_initial = x(:);
            [obj.x_normalized, obj.xmin,obj.xmax]  = AbstractPolynomial.normalize_x(x(:));
        end
        function fill_vander_matrix(obj,x,degree)
            fill_x(obj,x);
            [obj.vandermatrix,obj.Q,obj.R] = eval_vandermatrix(obj,obj.x_normalized,degree);
        end
        function [V,Q,R] = eval_vandermatrix(obj,normalized_x,degree)
            % evaluates vandermatrix of polynomial of specified type
            V = zeros(numel(normalized_x), degree + 1);
            for i = 0:degree
                V(:, i + 1) = monomial(obj,degree,i,normalized_x);
            end
            if nargout>1
                [Q, R] = qr(V,"econ");
            end
        end
        function fill_derivative_vandermatrix(obj,x,degree)
            % fills object's vandermatrix as a derivative of the polynomial
            fill_x(obj,x);    
            [obj.vandermatrix,obj.R,obj.Q] = vandermatrix_derivative(obj,obj.x_normalized,degree);

        end
        function [V,Q,R] = vandermatrix_derivative(obj,x,degree)
            % evalautes derivative vandermatrix for the polynomial of specified degree
            n = degree - 1;
            V = zeros(numel(x), n + 1);
            for i = 0:n
                V(:, i + 1) = monomial_derivative(obj,degree,i + 1,x);
            end
            [Q, R] = qr(V,"econ");
        end
    end
    methods (Static)
        function [a,b] = scalers()
            a = AbstractPolynomial.left;
            b = AbstractPolynomial.right;
        end
        function [x_normalized,x_min,x_max] = normalize_x(x)
            [a,b] = AbstractPolynomial.scalers();
            if AbstractPolynomial.isnormalized(x)
                x_normalized = x;
                x_min = a;
                x_max = b;
                return
            end
            x_min = min(x); x_max = max(x);
            if ~issorted(x)
                x = sort(x);
            end
            x_normalized = a + (b-a)*(x - x_min)/(x_max - x_min);
        end
        function x = denormalize_x(normalized_x,x_min,x_max)
            [a,b] = AbstractPolynomial.scalers();
            x = x_min + (normalized_x - a)*(x_max - x_min)/(b-a);
        end
        function is_norm = isnormalized(x)
            [a,b] = AbstractPolynomial.scalers();
            is_norm = ~isempty(x) && issorted(x) && (x(1) == a) && (x(end) == b);
        end
    end
end