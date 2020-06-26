% performance: the strdist between an ideal seq and current seq at seqtime
% data is spike data:x:time y:neuron #
% seq_len: the length of the sequence, or the number of groups
% size_group: neuron numbers in each group
% start_neuron: the sequence start neuron number
% period:  the stimulation periods
% Aug 10,2016, Edited by Yina Wei
% Edited by Oscar Gonzalez

function [seq_time,performance]=compute_performance_spiketime2(data,seq_len,size_group,ideal_seq,start_neuron,period,start_time,end_time)

max_time_response=350; % Frame window size

seq_time=start_time:period:end_time; % time when the stimulation starts;

maxplot=min(length(seq_time),100);

for i=1:length(seq_time)
    tempdata=data(:,seq_time(i):seq_time(i)+max_time_response-1);   % select data for each testing/training windows
    time_group=[];
    spike_group=[];
    for j=1:seq_len
        indx=(start_neuron+size_group*(j-1)):(start_neuron+size_group*j-1);
        poptempdata=sum(tempdata(indx,:)); % population activity   
        winsize=50;
        kernel=gausswin(winsize);
        convdata=conv(poptempdata,kernel,'same');%conv(spikes(k,:),normpdf(1:10,0,1),'same');
        peaks=find(convdata==max(convdata)&convdata~=min(convdata));
        if length(peaks)~=0
            time_group=[time_group peaks(1)];
            spike_group=[spike_group j];
       end
    end
    [tmp,sorti]=sort(time_group);
    seq=spike_group(sorti);
    performance(i)=strdist(seq,ideal_seq);    
end

function sd=strdist(inp,ideal)
% Computes string distance between input sequence (inp) and ideal sequence (ideal)
  tmp=0;
  removei=[];
  
  % Remove from ideal that is not present in input
  for sii = 1:length(ideal)
    if (length(find(ideal(sii)==inp))==0)
      removei=[removei sii];
    end
  end
  ideal(removei)=[];
  
  sd=2*length(ideal);
  for sii = 1:length(inp)
      sd=sd-abs(sii-find(ideal==inp(sii)));
  end
