# PairTCR Web Interface

一个用于运行PairTCR流程的现代化Web界面，支持实时日志监控和后台任务执行。

## 功能特性

- 🌐 **直观的Web界面**: 基于Bootstrap 5的现代化响应式设计
- 🔄 **实时监控**: WebSocket实时日志输出和进度跟踪
- 🚀 **后台执行**: 任务在后台运行，无需保持浏览器窗口
- 📦 **一键下载**: 自动打包结果文件为ZIP格式
- 📱 **移动友好**: 支持手机和平板设备访问
- ⚙️ **参数配置**: 支持所有PairTCR流程参数的自定义设置

## 快速开始

### 1. 安装依赖

```bash
cd web
pip install -r requirements.txt
```

### 2. 启动服务器

```bash
python start_server.py
```

或者使用自定义端口：

```bash
python start_server.py --port 8080
```

### 3. 访问Web界面

打开浏览器访问：http://localhost:5000

## 使用方法

### 基本使用

1. **输入目录**: 输入包含FASTQ文件的目录路径
2. **设置参数**: 配置文件前缀、读取限制、线程数等参数
3. **启动流程**: 点击"Start Pipeline"开始分析
4. **监控进度**: 实时查看执行日志和进度条
5. **下载结果**: 流程完成后一键下载压缩包

### 高级选项

- **MiXCR JAR路径**: 自定义MiXCR软件位置
- **C版本加速**: 使用编译的C程序提高性能
- **强制重启**: 删除已有结果重新开始

### 参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| Input Directory | FASTQ文件目录 | 必需 |
| File Prefix | 文件前缀 | TCR_TSO_18 |
| Read Limit | 最大读取数量 | 100,000 |
| Threads | CPU线程数 | 4 |
| MiXCR JAR | MiXCR程序路径 | ../scripts/mixcr.jar |
| Use C Version | 使用C版本加速 | False |
| Force Restart | 强制重新开始 | False |

## 目录结构

```
web/
├── app.py              # 主Flask应用
├── pipeline_runner.py  # Web版本的流程运行器
├── start_server.py     # 服务器启动脚本
├── requirements.txt    # Python依赖
├── README.md           # 本文档
├── templates/          # HTML模板
│   ├── index.html      # 主页面
│   └── job.html        # 任务监控页面
├── static/             # 静态文件（自动创建）
└── results/            # 结果文件（自动创建）
```

## 技术架构

- **后端**: Flask + Flask-SocketIO
- **前端**: Bootstrap 5 + Vanilla JavaScript
- **实时通信**: WebSocket (Socket.IO)
- **任务管理**: Python threading
- **文件处理**: 自动ZIP压缩

## 命令行选项

启动服务器时可以使用以下选项：

```bash
python start_server.py [选项]

选项:
  --host HOST       绑定主机地址 (默认: 0.0.0.0)
  --port PORT       端口号 (默认: 5000)
  --debug           调试模式
  --check-only      仅检查依赖，不启动服务器
```

## 生产部署

### 使用Gunicorn

```bash
pip install gunicorn
gunicorn --worker-class eventlet -w 1 --bind 0.0.0.0:5000 app:app
```

### 使用Nginx反向代理

```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://127.0.0.1:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }

    location /socket.io/ {
        proxy_pass http://127.0.0.1:5000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }
}
```

## 故障排除

### 常见问题

1. **端口被占用**
   ```bash
   python start_server.py --port 8080
   ```

2. **缺少依赖包**
   ```bash
   pip install -r requirements.txt
   ```

3. **无法找到MiXCR**
   - 检查MiXCR JAR文件路径
   - 确保文件权限正确

4. **C版本无法使用**
   - 确保已编译C程序
   - 检查可执行文件权限

### 日志调试

启用调试模式查看详细日志：

```bash
python start_server.py --debug
```

## 安全注意事项

⚠️ **重要**: 此Web界面主要用于内网环境，如需公网部署请注意：

1. 添加用户认证机制
2. 限制文件访问权限  
3. 使用HTTPS协议
4. 配置防火墙规则

## 许可证

本项目继承PairTCR的许可证条款。 