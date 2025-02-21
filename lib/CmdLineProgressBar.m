classdef CmdLineProgressBar < handle
    % CMDLINEPROGRESSBAR A class for command-line progress bar notification.
    %
    %   This class provides a simple way to display a progress bar in the
    %   command line during iterative operations.
    %
    % Example usage:
    %   pb = CmdLineProgressBar('Processing: ');
    %   for k = 1:10
    %       pb.print(k, 10);
    %       % Perform tasks here
    %       pause(0.5); % Simulating task delay
    %   end
    %
    % Project: statFEM-Recon
    % Author: Vahab Narouie
    % License: GNU GPL v3.0 (see LICENSE file for details)
    % Original inspiration from Itamar Katz, itakatz@gmail.com
    
    properties
        last_msg_len = 0; % Length of the last progress message displayed
    end
    
    methods
        % Constructor
        function obj = CmdLineProgressBar(msg)
            % CMDLINEPROGRESSBAR Constructor to initialize the progress bar.
            %
            % Inputs:
            %   msg - A message string displayed before the progress bar.
            %
            % Example:
            %   pb = CmdLineProgressBar('Loading: ');
            
            fprintf('%s', msg); % Print the initial message
        end
        
        % Print method
        function print(obj, n, tot)
          % PRINT Updates the progress bar with the current progress.
          %
          % Inputs:
          %   n   - Current iteration or progress count.
          %   tot - Total number of iterations or tasks.
          %
          % Example:
          %   pb.print(3, 10);
          
          % Erase the last progress message
          fprintf('%s', char(8 * ones(1, obj.last_msg_len)));
          
          % Create and display the new progress message
          info_str = sprintf('%d/%d', n, tot);
          fprintf('%s', info_str);
          
          % If this is the last step, print a newline for clean formatting
          if n == tot
            fprintf('\n');
          end
          
          % Update the length of the last message
          obj.last_msg_len = length(info_str);
        end
        
        % Destructor
        function delete(obj)
          % DELETE Cleans up the progress bar when the object is deleted.
          fprintf('\n'); % Ensure a newline at the end
        end
        
    end
    
end
