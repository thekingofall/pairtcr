# PairTCR Web Interface 使用指南

## 快速开始

### 1. 安装Web依赖

```bash
cd PairTCR/web
pip install -r requirements.txt
```

### 2. 启动Web服务器

有两种方式启动：

**方式一：从web目录启动**
```bash
cd web
python start_server.py
```

**方式二：从项目根目录启动（推荐）**
```bash
cd PairTCR
python start_web.py
```

### 3. 访问Web界面

浏览器打开：http://localhost:5000

## 使用步骤

### 第一步：提交任务
1. 在"Input Directory"中输入FASTQ文件所在目录的完整路径
2. （可选）设置"Output Directory"指定结果输出路径，如果留空则会在输入目录同级创建以文件前缀命名的文件夹
3. 可选择修改其他参数（文件前缀、读取限制、线程数等）
4. 点击"Advanced Options"可配置高级选项
5. 点击"Start Pipeline"提交任务

### 第二步：管理后台任务
- 任务提交后会显示唯一的任务ID，请妥善保存
- 可以选择"立即查看"跳转到监控页面，或"继续提交任务"留在主页面
- 在主页面的"查看任务"区域，输入任务ID可随时查看任务状态
- 任务会在后台持续运行，可以关闭浏览器

### 第三步：监控进度
- 实时查看执行日志和进度条
- 支持自动滚动和暂停功能
- 显示详细的任务参数和运行状态

### 第四步：下载结果
- 任务完成后，监控页面会显示"Download Results"按钮
- 点击即可下载包含所有结果文件的ZIP压缩包

## 参数说明

| 参数 | 必需 | 默认值 | 说明 |
|------|------|--------|------|
| Input Directory | ✓ | - | 包含FASTQ文件的目录路径 |
| Output Directory | | 输入目录同级/前缀名 | 结果输出目录路径 |
| File Prefix | | TCR_TSO_18 | 输出文件前缀 |
| Read Limit | | 100,000 | 最大处理reads数量 |
| Threads | | 4 | MiXCR使用的CPU线程数 |
| MiXCR JAR | | ../scripts/mixcr.jar | MiXCR程序路径 |
| Use C Version | | 否 | 使用C版本预处理（更快） |
| Force Restart | | 否 | 强制重新开始（删除现有结果） |

## 常见问题

**Q: 提示"Input directory does not exist"**
A: 请确保输入的路径存在且可访问，使用绝对路径更可靠

**Q: 页面显示"缺少依赖包"**  
A: 运行 `pip install -r requirements.txt` 安装所需依赖

**Q: 端口5000被占用**
A: 使用其他端口：`python start_server.py --port 8080`

**Q: 任务运行很慢**
A: 可以尝试启用"Use C Version"选项，或增加线程数

**Q: 无法下载结果**
A: 检查磁盘空间是否充足，确保有写入权限

**Q: 忘记了任务ID怎么办**
A: 在浏览器的历史记录中查找，或查看Recent Jobs列表

**Q: 任务ID查询显示"未找到"**
A: 请检查任务ID是否输入正确，确保是完整的UUID格式

**Q: 任务提交后想继续提交其他任务**
A: 点击成功提示框中的"继续提交任务"按钮，任务ID会自动填入查看框

## 目录结构

运行后web目录下会创建：
```
web/
├── results/          # 任务结果目录
│   ├── job-uuid-1/   # 各个任务的输出
│   ├── job-uuid-2/
│   └── *.zip         # 下载用的压缩包
└── static/           # 静态文件目录
```

## 生产部署建议

对于多用户或长期使用，建议：

1. **使用Gunicorn运行**
   ```bash
   pip install gunicorn
   gunicorn --worker-class eventlet -w 1 --bind 0.0.0.0:5000 app:app
   ```

2. **配置Nginx反向代理**
3. **定期清理results目录中的旧文件**
4. **考虑添加用户认证机制**

## 技术支持

如有问题，请检查：
1. Python版本是否为3.7+
2. 所有依赖是否正确安装
3. MiXCR文件是否存在且可执行
4. 输入数据格式是否正确

详细技术文档请参考：[README.md](README.md) 