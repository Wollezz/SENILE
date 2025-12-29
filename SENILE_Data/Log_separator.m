clear; clc;
folder = "Roccaraso2025";
addpath("folder");
fileList = dir(fullfile(folder, '*.csv'));

% Parametri
windowSize = 3000;  
maxRange = 15;       

for i = 1:length(fileList)
    fname = fileList(i).name;
    fullfilepath = fullfile(folder, fname);
    
    % Caricamento dati
    T = readtable(fullfilepath);
    dataMatrix = table2array(T);
    
    % Selezione asse accelerazione
    if startsWith(fname, "VN100")
        accy = dataMatrix(:, 12); 
    elseif startsWith(fname, "LSM6D")
        accy = dataMatrix(:, 3);
    else
        continue;
    end

    %% 1. Individuazione Ignition (Cut iniziale)
    % Nota: hai impostato < -20 (immagino per orientamento sensore)
    ignition_idx = find(accy < -20, 1, 'first');
    if isempty(ignition_idx), continue; end
    
    pre_cut_idx = max(1, ignition_idx - 100);

    %% 2. Individuazione Parte Post-Lancio (Ricerca segmento più lungo)
    accy_post = accy(ignition_idx:end);
    
    movingMax = movmax(accy_post, [windowSize-1, 0]);
    movingMin = movmin(accy_post, [windowSize-1, 0]);
    movingRange = movingMax - movingMin;
    
    % Vettore booleano: 1 se il punto fa parte di una finestra stabile
    isStable = movingRange <= maxRange;
    
    % --- LOGICA PER TROVARE IL SEGMENTO PIÙ LUNGO ---
    % Identifichiamo i blocchi di '1' consecutivi in isStable
    d = diff([0; isStable; 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    
    if ~isempty(starts)
        % Calcoliamo la durata di ogni segmento trovato
        durations = ends - starts + 1;
        
        % Troviamo l'indice del segmento più lungo
        [maxDuration, maxIdx] = max(durations);
        
        % Indici relativi alla sezione post-lancio
        startRel = starts(maxIdx) - windowSize + 1; % Torniamo all'inizio della finestra
        endRel = ends(maxIdx);
        
        % Controllo sicurezza indici
        startRel = max(1, startRel);
        
        % --- RICONVERSIONE INDICI PER IL PLOT TOTALE ---
        globalStartStable = ignition_idx + startRel - 1;
        globalEndStable = ignition_idx + endRel - 1;
        
        %% 3. Plot Completo
        figure('Name', ['Analisi: ' fname], 'NumberTitle', 'off');
        hold on;
        
        % Fase di Volo completa (Sfondo grigio)
        plot(1:length(accy), accy, 'Color', [0.8 0.8 0.8], 'DisplayName', 'Segnale Integrale');
        
        % Parte tagliata all'inizio (Blu)
        plot(1:pre_cut_idx, accy(1:pre_cut_idx), 'b', 'LineWidth', 1, 'DisplayName', 'Tagliato (Pre-Lancio)');
        
        % La parte stabile più lunga individuata (Rosso)
        plot(globalStartStable:globalEndStable, accy(globalStartStable:globalEndStable), 'r', 'LineWidth', 2, 'DisplayName', 'Parte Stabile Più Lunga');
        
        xline(ignition_idx, '--k', 'Ignition', 'LabelVerticalAlignment', 'bottom');
        
        title(sprintf('File: %s\nSegmento stabile più lungo: %d samples', fname, maxDuration), 'Interpreter', 'none');
        ylabel('Accelerazione [m/s^2]');
        xlabel('Samples');
        legend('Location', 'best');
        grid on;
        hold off;

        
        % Export 
        outputFolder = folder + "_separati"; 
        
        % Crea la cartella se non esiste già
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end
        
        fname_no_ext = fname(1:end-4); % Rimuove .csv dal nome originale
        
        writetable(T(1:pre_cut_idx, :), fullfile(outputFolder, [fname_no_ext '_pre.csv']));
        writetable(T(globalStartStable:globalEndStable, :), fullfile(outputFolder, [fname_no_ext '_post.csv']));
        
    else
        fprintf('File %s: Nessuna parte stabile trovata.\n', fname);
        figure('Name', ['Analisi Fallita: ' fname]);
        plot(1:length(accy), accy, 'Color', [0.8 0.8 0.8]); hold on;
        plot(1:pre_cut_idx, accy(1:pre_cut_idx), 'b');
        xline(ignition_idx, '--k', 'Ignition');
        title(['Nessuna zona stabile trovata in: ', fname], 'Interpreter', 'none');
        grid on;
    end
end