---
title: 博客维护与写作备忘录
date: 2026-01-28 16:31:25
tags: [备忘录, Hexo]
categories: [博客日志]
---

这篇文档用于记录本博客的维护流程、写作规范以及常用配置，防止遗忘。

<!-- more -->

## 🛠️ 日常操作工作流

### 1. 新建文章

在 VS Code 终端执行：

```bash
hexo new "文章标题"
```

或者直接在 `source/_posts/` 下新建 `.md` 文件。

### 2. 本地预览

修改完文章或配置后，执行“三连击”清理缓存并预览：

```bash
hexo clean
hexo g
hexo s
```

访问：[http://localhost:4000](http://localhost:4000)

---

## 📝 Markdown 常用语法速查

### 基础格式

- **加粗**: `**文字**` -> **文字**
- *斜体*: `*文字*` -> *文字*
- 引用: `> 这是一段引用`
- 链接: `[显示文字](链接地址)`

### 插入图片

1. 把图片放在 `source/img/` 文件夹下（例如 `test.jpg`）。
2. 在文章中引用：

```markdown
![图片描述](/img/test.jpg)
```

*(注意：路径最前面有个斜杠 `/`)*

### 代码块

使用三个反引号包裹：

```python
print("Hello World")
```

### 标签插件

1. 彩色提示块

```markdown
{% note [颜色] [no-icon] %}
这里写你的内容...
支持 **Markdown** 和 $ E=mc^2 $
{% endnote %}
```

效果：

{% note [颜色] [no-icon] %}
这里写你的内容...
支持 **Markdown** 和 $ E=mc^2 $
{% endnote %}

- **颜色选项：default(灰色)，primary(紫色)，info(蓝色)，success (绿色)，warning (黄色)，danger (红色)**

1. 标签页盒子

```markdown
{% tabs test-id %}
<!-- tab 这里的标题1 -->
这里是第一个标签页的内容。
比如放 Python 代码。
<!-- endtab -->

<!-- tab 这里的标题2 -->
这里是第二个标签页的内容。
比如放 C++ 代码。
<!-- endtab -->

<!-- tab 这里的标题3 -->
这里是第三个内容。
<!-- endtab -->
{% endtabs %}
```

效果：

{% tabs test-id %}
<!-- tab 这里的标题1 -->
这里是第一个标签页的内容。
比如放 Python 代码。
<!-- endtab -->

<!-- tab 这里的标题2 -->
这里是第二个标签页的内容。
比如放 C++ 代码。
<!-- endtab -->

<!-- tab 这里的标题3 -->
这里是第三个内容。
<!-- endtab -->
{% endtabs %}

---

## 📐 数学公式 (LaTeX)

博客已配置 MathJax 渲染，支持行内和块级公式。

### 行内公式

写法：`$ E=mc^2 $`
效果：$ E=mc^2 $

### 块级公式

写法：

```latex
$$
\begin{equation*}
  \sum_{i=1}^n i = \frac{n(n+1)}{2}
\end{equation*}
$$
```

预览：

$$
\begin{equation*}
  \sum_{i=1}^n i = \frac{n(n+1)}{2}
\end{equation*}
$$

### 多行对齐 (Align)

写法：

```latex
$$
\begin{align}
  f(x) &= (x+a)(x+b) \\
       &= x^2 + (a+b)x + ab
\end{align}
$$
```

预览：
$$
\begin{align}
  f(x) &= (x+a)(x+b) \\
       &= x^2 + (a+b)x + ab
\end{align}
$$

{% note danger %}
**经过测试发现使用行内公式时内容必须和美元符号间有一个空格，行间公式不能直接使用 `$$...$$` 的形式，如果不按此规范会导致内容无法正常渲染。**
{% endnote %}

---

## 🏷️ 文章头部 (Front-matter) 模板

每次新建文章，头部必须包含以下信息：

```yaml
---
title: 文章标题
date: 2024-01-28 14:30:00
updated: 2024-01-28 14:30:00  # 可选，更新时间
tags: 
  - 标签1
  - 标签2
categories: 
  - [分类A, 子分类B]
top_img: /img/banner.jpg      # 可选，文章顶部背景图
cover: /img/cover.jpg         # 可选，首页显示的缩略图
---
```
