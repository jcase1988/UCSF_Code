function [banks, bank_labels] = default_banks(subj)

if strcmp(subj,'EC70')
    banks = [(1:16:256)' (16:16:256)'];
    banks = [banks ;
        257 262;  %OFA (6)
        263 268;  %OFP (6)
        269 272;  %TP  (4)
        273 276;  %AST (4)
        277 280;  %MST (4)
        281 284;  %PST (4)
        285 288;  %HDA (4)
        289 292;  %HDP (4)
        293 296]; %CD  (4)

elseif strcmp(subj,'EC77')
        banks = [(1:16:256)' (16:16:256)'];
        banks = [banks ;
            257 320];
                
    %setup labels for the beginning of each strip/depth
        bank_labels = {};
        bank_labels(1:256) = {'Grid'};
        bank_labels(257:260) = {'AOF'};
        bank_labels(261:264) = {'POF'};
        bank_labels(265:270) = {'ITG'};
        bank_labels(271:276) = {'Pos-Temp'};
        bank_labels(277:280) = {'AST'};
        bank_labels(281:284) = {'MST'};
        bank_labels(285:288) = {'PST'};
        bank_labels(289:292) = {'Insula'};
        bank_labels(293:296) = {'sACC'};
        bank_labels(297:300) = {'iACC'};
        bank_labels(301:304) = {'Amyg'};
        bank_labels(305:308) = {'Hippo'};
        
elseif strcmp(subj,'EC82')
        banks = [(1:16:256)' (16:16:256)'];
        banks = [banks ;
                257 276 ; 
                277 314];
            
        bank_labels = {};
        bank_labels(1:256) = {'Grid'};
        bank_labels(257:276) = {'FG'};
        bank_labels(277:280) = {'Amyg'};
        bank_labels(281:284) = {'Hippo'};
        bank_labels(285:288) = {'AST'};
        bank_labels(289:292) = {'MST'};
        bank_labels(293:296) = {'PST'};
        bank_labels(297:302) = {'ITG'};
        bank_labels(303:306) = {'F1'};
        bank_labels(307:310) = {'sOFC'};
        bank_labels(311:314) = {'iOFC'};
    
elseif strcmp(subj,'EC71')
        banks = [(1:16:128)' (16:16:128)'];
        banks = [banks ;
                129 164];
        
        bank_labels = {};
        bank_labels(1:64) = {'OFC'};
        bank_labels(65:128) = {'Grid'};
        bank_labels(129:132) = {'AST'};
        bank_labels(133:136) = {'MST'};
        bank_labels(137:140) = {'PST'};
        bank_labels(141:144) = {'AD'};
        bank_labels(145:148) = {'HD'};
        bank_labels(149:152) = {'AID'};
        bank_labels(153:156) = {'PID'};
        bank_labels(157:160) = {'SCD'};
        bank_labels(161:164) = {'ICD'};
        
        
end

end