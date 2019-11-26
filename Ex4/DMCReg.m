classdef DMCReg < handle
    %DMCReg Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LAMBDA
        M
        Mp
        K
        psi
        N
        Gz
        D
        Nu
        u_prev
        yzad
        deltaup
        deltazp
        Mpz
        z_prev
        Dz
    end
    
    methods
        function obj = DMCReg(Gz, D, N, Nu, lambda)
            %DM Construct an instance of this class
            %   Detailed explanation goes here
            obj.N = N;
            obj.Gz = Gz;
            obj.D = D;
            obj.Nu = Nu;
            s = step(Gz, D/Gz.Ts);
            s = s(:,:,1);
            obj.psi = eye(N);
            obj.LAMBDA = eye(Nu)*lambda;


            % Wyznaczanie macierzy M
            obj.M = zeros(N, Nu);
            for i=1:1:Nu
                obj.M(i:N,i)=s(1:N-i+1)';
            end


            % Wyznaczanie macierzy Mp
            obj.Mp= zeros(N, D-1);
            for i=1:N
               for j=1:D-1
                  if i+j<=D
                     obj.Mp(i,j)=s(i+j)-s(j);
                  else
                     obj.Mp(i,j)=s(D)-s(j);
                  end   
               end
            end

            % Wektor wzmocnieñ - wyznaczany raz (offline)
            obj.K=(obj.M'*obj.psi*obj.M+obj.LAMBDA)^(-1)*obj.M'*obj.psi;
        end
        
        function [] = reset(obj,u_p)
            obj.u_prev = u_p;
            obj.deltaup=zeros(1,obj.D-1)';
        end
        
        function [] = setValue(obj, yzad)
            yzad(1:obj.N)= yzad;
            obj.yzad=yzad';
        end
        
        function u = countValue(obj,y)
            %count actual value of output.
            %   Detailed explanation goes here
            
            % aktualizacja wektora aktualnej wartoœci wyjœcia
            yk=ones(obj.N,1)*y;

            % wyliczenie nowego wektora odpowiedzi swobodnej
            y0=yk+obj.Mp*obj.deltaup;

            % wyliczenie wektora zmian sterowania
            deltauk=obj.K*(obj.yzad-y0);

            % prawo regulacji
            u=obj.u_prev+deltauk(1);
            obj.u_prev = u;

            % aktualizacja poprzednich zmian sterowania
            obj.deltaup=[deltauk(1) obj.deltaup(1:end-1)']';
        end
    end
end

