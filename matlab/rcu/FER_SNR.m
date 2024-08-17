close all; clear; clc;

SNR_dB = [2.5, 3, 3.5, 4];

k  = 64;
n  = 128;
pX = [0.5 0.5];
X  = [-1 1];

format long;
rcu = zeros(1,length(SNR_dB));

for i_snr = 1:length(SNR_dB)
    rcu_result = RCU_warpper(k, n, pX, X, SNR_dB(i_snr));
    rcu(i_snr) = rcu_result(4);
end

dld_bam_ls_1e4 = [0.108342, 0.0300661, 0.00469395, 0.00081];
dld_bam_ls_2e4 = [0.0834028, 0.019658, 0.0035624, 0.00042];
figure();

semilogy(SNR_dB, dld_bam_ls_1e4, '-o', 'DisplayName', "BAM list size = 1e4", 'LineWidth', 1.5);
hold on;
semilogy(SNR_dB, dld_bam_ls_2e4, '-o', 'Color', 'r', 'DisplayName', "BAM list size = 2e4", 'LineWidth', 1.5);
semilogy(SNR_dB, rcu, '--', 'Color', 'k', 'DisplayName', "RCU", 'LineWidth', 1.3);
hold off;
grid on;
xlabel('SNR (dB)', 'Interpreter', 'latex', 'FontSize', 17);
ylabel('Frame Error Rate', 'Interpreter', 'latex', 'FontSize', 17);

ylim([1e-7, 1e-1]);
legend('Interpreter', 'latex', 'FontSize', 15);