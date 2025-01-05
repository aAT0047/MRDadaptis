import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd

# 假设的输入和输出维度
n_features = 9 # 示例特征数量
n_targets = 8   # 示例目标数量

# 定义神经网络模型
class MultiTargetRegression(nn.Module):
    def __init__(self):
        super(MultiTargetRegression, self).__init__()
        self.fc1 = nn.Linear(n_features, 128)
        self.fc2 = nn.Linear(128, 64)
        self.fc3 = nn.Linear(64, n_targets)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)
        return x

# 实例化模型、损失函数和优化器
model = MultiTargetRegression()
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# 加载元数据
# torch.manual_seed(42)  # 为了结果可复现
# inputs = torch.rand(100, n_features)  # 100个样本，每个样本n_features个特征
# targets = torch.rand(100, n_targets)  # 100个样本的目标值
# 加载元数据
df = pd.read_csv(r"C:\Users\1\Desktop\code_ppg\Autosvp\finaF1.csv")
df = df.iloc[:, 1:-3]
df = df.dropna()

# print(df)
# 分离特征和标签
feature = df.iloc[:, 2:-8]  # 所有行，除了最后8列
target = df.iloc[:, -8:]    # 所有行，只有最后8列

# 将NumPy数组转换为PyTorch张量
inputs = torch.tensor(feature.values, dtype=torch.float32)
# print(inputs)
targets = torch.tensor(target.values, dtype=torch.float32)
# 创建数据加载器
dataset = TensorDataset(inputs, targets)
train_loader = DataLoader(dataset, batch_size=4, shuffle=True)

# 训练模型
epochs =31000
model.train() 
for epoch in range(epochs):
    for batch_idx, (data, target) in enumerate(train_loader):
        # 前向传播
        outputs = model(data)
        loss = criterion(outputs, target)

        # 反向传播和优化
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if batch_idx % 10 == 0:
            print(f'Epoch [{epoch+1}/{epochs}], Step [{batch_idx+1}/{len(train_loader)}], Loss: {loss.item():.4f}')

# 模型评估（这里我们简单地在训练集上评估）
model.eval()
with torch.no_grad():
    predictions = model(inputs)
    loss = criterion(predictions, targets)
    print(f'\nEvaluation Loss: {loss.item():.4f}')
# # 保存模型的状态字典
# torch.save(model.state_dict(), 'multi_target_regression_model.pth')
import os

# 打印当前工作目录
current_directory = os.getcwd()
print("当前工作目录:", current_directory)

torch.save(model.state_dict(), r'C:\Users\1\Desktop\code_ppg\Autosvp\multi_target_regression_model1.pth')
