% load transmission line frequency-dependent RLGC parameters
load stripline_model

% transmission line length (in metres)
d = 0.8;

% create transmission line model
tline = rlgc(R,L,G,C,d,f);

% bondwire package
len = 1e-3;
pad1 = admittance(j*2*pi*f*0.45e-12);
bwire = impedance(j*2*pi*f*len*9.6e-7);
pad2 = admittance(j*2*pi*f*0.45e-12);
pkg = series(pad1,series(bwire,pad2));

% source impedance
source = impedance(120*ones(1,length(f)));

% termination
termination = admittance(ones(1,length(f))./80);

% channel is the series connection of the source, package, transmission
% line, package, and termination
channel = series(source,series(series(pkg,series(tline,pkg)),termination));

% frequency domain response
Hchannel = 1./channel.A;

% time domain response
[h,t,hstep] = freq2impulse(Hchannel,f);

hstep = hstep(500:2:end);
h = h(500:2:end);
t = t(500:2:end);
Ts = t(2) - t(1);

% plot results in time and frequency domain
subplot(2,1,1);
plot(1e-9*f,20*log10(abs(Hchannel)),'linewidth',2);
set(gca,'xlim',[0 20],'fontsize',16);
grid on
xlabel('Frequency [GHz]');
ylabel('Mag Response (dB)');
subplot(2,1,2);
plot(t*1e9,hstep,'r-','linewidth',2);
set(gca,'xlim',[0 20],'ylim',[-.02 .42],'fontsize',16);
grid on
xlabel('t (ns)');
ylabel('Step Response')
