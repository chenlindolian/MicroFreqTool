function signal = simulateSignalFrom3DField(model, field, Params)

 omega = 2*pi*field;
 
 signal_model = model.model.*model.mask;
 
 N = length(Params.TEs);
 
 intra_axonal_index = (signal_model == 3);
 myelin_index = (signal_model == 2); 
 extra_axonal_index = (signal_model == 1);

 nb_pixel = sum(model.mask, 'all');

 AVF = sum(intra_axonal_index , 'all') / nb_pixel;
 MVF = sum(myelin_index , 'all') / nb_pixel;
 FVF = AVF+MVF;
 
for l = 1:N    
    intra_axonal = exp(1i*Params.TEs(l)*omega(intra_axonal_index));
    myelin = exp(1i*Params.TEs(l)*omega(myelin_index));
    extra_axonal = exp(1i*Params.TEs(l)*omega(extra_axonal_index));
    
    signal.intra_axonal(l) = sum(intra_axonal, 'all');
    signal.myelin(l) = sum(myelin, 'all');
    signal.extra_axonal(l) = sum(extra_axonal, 'all');
end

 signal.intra_axonal = AVF * Params.intra_axonal.proton_density * Params.intra_axonal.Mag * exp(-Params.TEs / Params.intra_axonal.T2).* signal.intra_axonal;
 signal.myelin = MVF * Params.myelin.proton_density * Params.myelin.Mag * exp(-Params.TEs / Params.myelin.T2).* signal.myelin;
 signal.extra_axonal = (1-FVF) * Params.extra_axonal.proton_density * Params.extra_axonal.Mag * exp(-Params.TEs / Params.extra_axonal.T2).*signal.extra_axonal;

 signal.total = signal.intra_axonal + signal.myelin + signal.extra_axonal;
 signal.total_normalized = signal.total ./ abs(signal.total(1));

end