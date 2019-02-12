function emailSimStatus(recipient,varargin)

%%%%% Email Simulation Status %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   emailSimStatus(RECIPIENT) sends the default message: "Your simulation
%   is finished" to RECIPIENT with a standard subject.
%
%   emailSimStatus(RECIPIENT,MESSAGE) allows you to send a custom MESSAGE
%
%   emailSimStatus(RECIPIENT,MESSAGE,SUBJECT) allows to to choose the
%   SUBJECT as well.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    subj = 'Simulation Status';
    mssg = 'Your simulation is finished';
elseif nargin == 2
    subj = 'Simulation Status';
    mssg = varargin{1};
elseif nargin == 3
    subj = varargin{2};
    mssg = varargin{1};
end

sender = 'jane.doe@gmail.com';
psswd  = 'password';

setpref('Internet','SMTP_Server','smtp.gmail.com'); 
setpref('Internet','E_mail',sender);
setpref('Internet','SMTP_Username',sender);
setpref('Internet','SMTP_Password',psswd);
 
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(recipient,subj,mssg)

