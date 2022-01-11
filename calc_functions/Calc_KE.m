function KE = Calc_KE(Q_save,LastDomainZindex)

% Compute GW kinetic energy

% INPUTS: 
% Q_save -  4-D array (inclusive of ghost cells/ sponge)
% LastDomainZindex - last (max) vertical index of the computational domain.
            

KE = squeeze(0.5*(Q_save(3:LastDomainZindex,3:end-2,2,:).^2+Q_save(3:LastDomainZindex,3:end-2,3,:).^2)./Q_save(3:LastDomainZindex,3:end-2,1,:));

end