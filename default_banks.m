function banks = default_banks()

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
end