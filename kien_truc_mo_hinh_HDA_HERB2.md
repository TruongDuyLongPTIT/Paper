# Kiến trúc mô hình dự đoán Herb-Disease Association trên HERB 2.0
## Backbone BG-HGNN (có trọng số bằng chứng) + Nhánh Cross-channel Attention (kiểu iCAM-Net)

> Tài liệu này giải thích từng bước tính toán một cách chi tiết nhất có thể — mỗi công thức đều đi kèm (1) giải thích bằng lời trước khi đưa công thức, (2) chính công thức, (3) một ví dụ số cụ thể để bạn tự tay dò lại từng bước. Mục tiêu là đọc xong có thể tự giải thích lại cho người khác mà không cần nhìn công thức gốc.

---

## 0. Bối cảnh, mục tiêu và lý do chọn kiến trúc này

### 0.1. Bài toán đang giải quyết

Đầu vào là một **đồ thị tri thức** (knowledge graph) tên HERB 2.0, mô tả thế giới của thuốc thảo dược (herb), bệnh (disease), và các khái niệm liên quan. Đồ thị này có:

- **9 loại thực thể** (node type): herb, disease, ingredient (thành phần hoạt chất), target protein, symptom (triệu chứng), syndrome (hội chứng), meridian (kinh mạch), prescription (bài thuốc), GO term/pathway (con đường sinh học).
- **28 loại quan hệ** (edge type): ví dụ herb–ingredient (herb chứa thành phần nào), disease–target (bệnh liên quan protein nào), herb–meridian (herb đi vào kinh mạch nào), herb–symptom, v.v.
- Mỗi cạnh cụ thể trong đồ thị còn có một **trọng số bằng chứng liên tục** — con số phản ánh độ tin cậy của quan hệ đó, chia thành 4 tầng: thử nghiệm lâm sàng (đáng tin nhất) > phân tích tổng hợp (meta-analysis) > thí nghiệm in-vitro/in-vivo > suy luận thống kê từ máy tính (kém tin cậy nhất).

**Mục tiêu**: cho một cặp (herb, disease) chưa từng được ghi nhận có liên hệ, dự đoán xác suất chúng thực sự liên quan tới nhau — và quan trọng hơn, giải thích được **vì sao** mô hình nghĩ vậy (thành phần hoạt chất nào của herb, protein nào của disease là "thủ phạm chính" đứng sau dự đoán).

### 0.2. Vì sao không dùng thẳng kiến trúc có sẵn (HGHDA, iCAM-Net, R-GCN thuần)

- **HGHDA và iCAM-Net**: cứng hóa kiến trúc quanh đúng 4 loại thực thể (herb, ingredient, protein, disease) và 2 quan hệ chính (herb-ingredient, disease-protein). Không có chỗ để nhét thêm symptom, meridian, prescription, pathway — tức là bỏ phí phần lớn thông tin phong phú của HERB 2.0, và không có khái niệm trọng số bằng chứng liên tục (chỉ dùng ma trận nhị phân 0/1).
- **R-GCN gốc (Schlichtkrull et al.)**: xử lý được đúng tình huống nhiều loại thực thể/nhiều loại quan hệ, nhưng với 28 loại quan hệ, số tham số tăng tuyến tính theo số quan hệ (mỗi quan hệ một ma trận trọng số riêng W_r) — nguy cơ quá khớp (overfit) rất cao với các quan hệ có ít cạnh dữ liệu.
- **BG-HGNN (Su et al. 2024)**: giải quyết đúng vấn đề "nhiều loại quan hệ" của R-GCN bằng cách không cần một ma trận trọng số riêng cho từng quan hệ nữa — số tham số không phụ thuộc vào số lượng quan hệ. Đây là lý do chọn nó làm **backbone** (xương sống) xử lý toàn bộ 9 thực thể/28 quan hệ.
- Nhưng BG-HGNN gộp mọi thông tin vào một không gian chung, nên **mất khả năng tách riêng "component nào tương tác mạnh với protein nào"** cho một cặp herb-disease cụ thể — đây chính là điểm mạnh khiến iCAM-Net diễn giải được ở mức phân tử (case study Angelica sinensis + 5′-GMP–APOE). Vì vậy cần thêm **nhánh attention riêng** chạy song song, chỉ tập trung vào đúng 2 quan hệ herb-ingredient/disease-target, để không đánh mất khả năng diễn giải này.

→ Kết luận thiết kế: **BG-HGNN làm xương sống xử lý toàn bộ đồ thị (có trọng số bằng chứng) + nhánh attention riêng kiểu iCAM-Net để diễn giải phân tử**, hai nhánh chạy song song rồi gộp lại.

---

## 1. Sơ đồ tổng thể luồng dữ liệu

```
                         ┌─────────────────────────────┐
                         │   KG HERB 2.0 (đầu vào)      │
                         │ 9 loại thực thể, 28 quan hệ, │
                         │ trọng số bằng chứng w_uv     │
                         └──────────────┬──────────────┘
                                        │
                 ┌──────────────────────┴───────────────────────┐
                 │                                               │
                 ▼                                               ▼
   ┌─────────────────────────────┐            ┌───────────────────────────────────┐
   │   NHÁNH A — BG-HGNN         │            │  NHÁNH B — Cross-channel Attention │
   │   backbone toàn KG          │            │  (chỉ herb-ingredient/            │
   │   (mọi thực thể, mọi quan   │            │   disease-target)                 │
   │   hệ, có trọng số w_uv)     │            │  (Q/K/V hai chiều, kiểu iCAM-Net)  │
   └──────────────┬───────────────┘            └────────────────┬───────────────────┘
                  │ E_h^BG, E_d^BG                                │ E_h^attn, E_d^attn
                  │ (+ embedding mọi node khác)                   │ (+ attention map A_cp, A_pc)
                  └──────────────────┬─────────────────────────────┘
                                     ▼
                    ┌────────────────────────────────────┐
                    │   NHÁNH C — Fusion Module           │
                    │   combine(E^BG, E^attn) → E_final   │
                    └──────────────────┬───────────────────┘
                                       ▼
                    ┌────────────────────────────────────────────┐
                    │   NHÁNH D — Multi-task heads                │
                    │   Task chính: HDA score = σ(E_h · E_d^T)     │
                    │   Task phụ: Ingredient-Target Association    │
                    │   (dùng F_C, F_P từ nhánh B)                  │
                    └────────────────────────────────────────────┘
```

Bốn khối A, B, C, D được giải thích chi tiết ở mục 2, 3, 4, 5 dưới đây, mỗi bước có ví dụ số minh họa.

---

## 2. NHÁNH A — BG-HGNN backbone có trọng số bằng chứng

Nhánh này có 5 bước con, chạy tuần tự: (2.1) chuẩn hóa feature, (2.2) mã hóa loại thực thể/quan hệ, (2.3) tính vector "bối cảnh quan hệ" có trọng số, (2.4) gộp 3 nguồn thông tin bằng Kronecker product + nén chiều, (2.5) lan truyền qua GNN có trọng số cạnh.

### 2.1. Chuẩn hóa đầu vào — Attribute Space Fusion

**Vấn đề cần giải quyết**: mỗi loại thực thể có feature gốc khác nhau về ý nghĩa và số chiều. Ví dụ herb có thể có feature là fingerprint hóa học tổng hợp từ các ingredient (giả sử 5 chiều), còn disease có feature là vector ngữ nghĩa từ MeSH (giả sử 3 chiều). Một GNN thông thường không thể xử lý cùng lúc 2 loại vector có số chiều khác nhau trong cùng một phép tính ma trận.

**Cách giải quyết**: tạo một "không gian chung" đủ rộng để chứa mọi loại kênh feature của mọi loại thực thể (gộp toàn bộ, không trùng lặp), sau đó với mỗi node, đặt giá trị thật vào đúng vị trí kênh của nó, còn các vị trí không thuộc về loại thực thể đó thì điền một giá trị mặc định (ví dụ 0, để mô hình biết "kênh này không áp dụng cho node này" — sự vắng mặt cũng là một tín hiệu).

**Công thức**:
```
X_τV = ⋃_τ X_τ                         (hợp mọi kênh feature của mọi loại thực thể)
x'_v = proj(x_v, X_τV)                  (điền giá trị đặc biệt cho kênh không tồn tại)
```

**Ví dụ số cụ thể**: giả sử toàn hệ thống chỉ có 2 loại thực thể — herb (feature gốc 5 chiều: f1..f5) và disease (feature gốc 3 chiều: g1..g3). Không gian chung sẽ có 8 chiều (5+3). Với herb "Đương quy" có feature gốc = [0.2, 0.5, 0.1, 0.8, 0.3], sau chiếu vào không gian chung sẽ thành:
```
x'_(Đương quy) = [0.2, 0.5, 0.1, 0.8, 0.3, 0, 0, 0]
                  └────── 5 kênh herb ──────┘ └ 3 kênh disease (không áp dụng) ┘
```
Với disease "Đái tháo đường" có feature gốc = [0.6, 0.4, 0.9], sau chiếu:
```
x'_(Đái tháo đường) = [0, 0, 0, 0, 0, 0.6, 0.4, 0.9]
                       └── 5 kênh herb (không áp dụng) ──┘ └ 3 kênh disease ┘
```
Giờ cả hai vector cùng nằm trong không gian 8 chiều, có thể đưa vào cùng một phép tính.

### 2.2. Mã hóa loại thực thể và loại quan hệ — Type Encoding

**Vấn đề cần giải quyết**: mô hình cần biết "node này thuộc loại nào" (herb, disease, ingredient...) và "cạnh này thuộc loại quan hệ nào" (herb-ingredient, disease-target...) — 9 loại thực thể và 28 loại quan hệ. Cách thông thường là dùng **one-hot encoding** (mỗi loại là một vector chỉ có đúng 1 số "1", còn lại toàn "0", ví dụ loại thứ 3 trong 9 loại là [0,0,1,0,0,0,0,0,0]).

**Vì sao không dùng one-hot**: theo thực nghiệm gốc của BG-HGNN, one-hot khiến mô hình hội tụ (converge) chậm hơn hẳn so với cách khác, vì one-hot tạo ra vector rất "thưa" (sparse — hầu hết là số 0), không giúp mô hình học được mối tương quan tinh tế giữa các loại.

**Cách giải quyết — dense random projection**: thay vì one-hot, gán cho mỗi loại thực thể/quan hệ một **vector ngẫu nhiên cố định** (random nhưng chọn một lần rồi giữ nguyên suốt quá trình train), lấy từ phân phối đều U(a,b). "Dense" nghĩa là mọi phần tử trong vector đều khác 0 (không thưa như one-hot).

**Công thức**:
```
node type encoding:  o_τv,i ∈ R^m,  mỗi phần tử rút ngẫu nhiên từ U(a,b)
edge type encoding:  o_τe,i ∈ R^n,  mỗi phần tử rút ngẫu nhiên từ U(a,b)
```

**Ví dụ số cụ thể**: giả sử m=4 (vector 4 chiều), U(a,b) = U(-1,1). Loại thực thể "herb" có thể được gán ngẫu nhiên một lần (rồi cố định mãi mãi):
```
o_herb       = [ 0.31, -0.72,  0.05,  0.88]
o_disease    = [-0.44,  0.19,  0.63, -0.27]
o_ingredient = [ 0.77,  0.02, -0.51,  0.40]
```
Ba vector này không có ý nghĩa "gần nhau hay xa nhau" như one-hot (nơi mọi loại đều cách đều nhau tuyệt đối) — chúng chỉ cần là các "chữ ký" (signature) đủ khác biệt để mô hình phân biệt được loại thực thể trong lúc học, và nhờ số chiều dày đặc, mô hình dễ tìm ra tương quan hơn.

### 2.3. Vector bối cảnh quan hệ có trọng số bằng chứng — điểm chèn quan trọng thứ nhất

Đây là bước quan trọng nhất cần hiểu kỹ vì đây là nơi trọng số bằng chứng được đưa vào lần đầu.

**Vấn đề cần giải quyết**: mỗi node v (ví dụ herb "Đương quy") có nhiều hàng xóm kết nối với nó qua nhiều loại quan hệ khác nhau (một số ingredient, một số disease, một số symptom...). Ta muốn tính ra **một vector duy nhất đại diện cho "bối cảnh quan hệ" của node v** — tức là "trung bình" các loại quan hệ mà v tham gia.

**Công thức gốc (trung bình đều, không phân biệt độ tin cậy)**:
```
ω_v = (1/|N(v)|) · Σ_{u∈N(v)} o_ψ(u,v)
```
Đọc thành lời: lấy tất cả hàng xóm u của v, với mỗi hàng xóm lấy vector loại-quan-hệ của cạnh (u,v) — gọi là o_ψ(u,v) — rồi cộng hết lại và chia đều cho số lượng hàng xóm. Đây đúng là phép "tính điểm trung bình cộng" bình thường.

**Vấn đề của công thức gốc**: mọi cạnh được coi bình đẳng như nhau. Một cạnh đến từ thử nghiệm lâm sàng (rất đáng tin) được cộng vào với trọng số y hệt một cạnh chỉ đến từ suy luận thống kê (kém tin hơn nhiều). Đây chính là điều ta muốn sửa.

**Công thức có trọng số (thay thế)**:
```
ω_v = ( Σ_{u∈N(v)} w_uv · o_ψ(u,v) ) / ( Σ_{u∈N(v)} w_uv + ε )
```
Đọc thành lời: thay vì cộng đều rồi chia cho số lượng, ta **nhân mỗi cạnh với trọng số bằng chứng của nó trước khi cộng**, rồi chia cho **tổng trọng số** (thay vì chia cho số lượng). Cạnh có trọng số cao đóng góp nhiều hơn vào tổng, cạnh có trọng số thấp đóng góp ít hơn. `ε` là một số rất nhỏ (ví dụ 1e-8) để tránh chia cho 0 nếu một node không có hàng xóm nào.

**Ví dụ số cụ thể — dò tay từng bước**: giả sử herb "Đương quy" có 3 hàng xóm, qua 3 loại quan hệ khác nhau:

| Hàng xóm u | Loại quan hệ | Vector loại quan hệ o_ψ(u,v) | Trọng số w_uv | Nguồn bằng chứng |
|---|---|---|---|---|
| Ferulic acid (ingredient) | herb-ingredient | [0.10, 0.60, -0.20, 0.30] | 1.00 | thử nghiệm lâm sàng |
| Đái tháo đường (disease) | herb-disease | [0.50, -0.10, 0.40, 0.05] | 0.85 | meta-analysis |
| Mất ngủ (symptom) | herb-symptom | [-0.30, 0.20, 0.60, -0.10] | 0.25 | suy luận thống kê |

**Bước 1 — trung bình đều (công thức gốc, để so sánh)**:
```
ω_v(gốc) = ( [0.10,0.60,-0.20,0.30] + [0.50,-0.10,0.40,0.05] + [-0.30,0.20,0.60,-0.10] ) / 3
         = [0.30, 0.70, 0.80, 0.25] / 3
         = [0.10, 0.233, 0.267, 0.083]
```
Chú ý: cạnh "Mất ngủ" (bằng chứng yếu nhất) vẫn đóng góp đúng 1/3 tỷ trọng — y hệt cạnh "Ferulic acid" (bằng chứng mạnh nhất).

**Bước 2 — trung bình có trọng số (công thức mới)**:

Tổng trọng số = 1.00 + 0.85 + 0.25 = 2.10

Nhân từng vector với trọng số tương ứng:
```
1.00 × [0.10, 0.60, -0.20, 0.30]  = [0.100, 0.600, -0.200, 0.300]
0.85 × [0.50, -0.10, 0.40, 0.05]  = [0.425, -0.085, 0.340, 0.0425]
0.25 × [-0.30, 0.20, 0.60, -0.10] = [-0.075, 0.050, 0.150, -0.025]
```
Cộng lại:
```
tổng = [0.100+0.425-0.075,  0.600-0.085+0.050,  -0.200+0.340+0.150,  0.300+0.0425-0.025]
     = [0.450, 0.565, 0.290, 0.3175]
```
Chia cho tổng trọng số (2.10):
```
ω_v(có trọng số) = [0.450/2.10, 0.565/2.10, 0.290/2.10, 0.3175/2.10]
                  = [0.214, 0.269, 0.138, 0.151]
```

**So sánh hai kết quả**:
```
ω_v(gốc, trung bình đều)    = [0.100, 0.233, 0.267, 0.083]
ω_v(có trọng số bằng chứng) = [0.214, 0.269, 0.138, 0.151]
```
Kết quả khác nhau rõ rệt — vector có trọng số nghiêng nhiều hơn về phía cạnh "Ferulic acid" (bằng chứng mạnh nhất, đóng góp gần một nửa tổng trọng số 1.00/2.10 ≈ 47.6%) so với cạnh "Mất ngủ" (chỉ đóng góp 0.25/2.10 ≈ 11.9%). Đây chính là hiệu ứng mong muốn: bằng chứng mạnh "có tiếng nói lớn hơn".

**Lưu ý quan trọng**: công thức này **không thêm tham số học được nào cả** — w_uv là số cố định tính sẵn trước khi train (từ 4 tầng bằng chứng), không phải trọng số mạng neural tự học. Điều này giữ nguyên lợi thế "số tham số không tăng theo số quan hệ" của BG-HGNN.

### 2.4. Gộp 3 nguồn thông tin — Kronecker product + Low-rank decomposition

**Vấn đề cần giải quyết**: sau các bước trên, mỗi node v có 3 vector riêng biệt: feature gốc x'_v (mục 2.1), vector loại thực thể o_τv,i (mục 2.2), và vector bối cảnh quan hệ có trọng số ω_v (mục 2.3, vừa tính ở trên). Cần gộp 3 vector này thành **một** vector duy nhất để đưa vào GNN. Cách đơn giản nhất là nối chúng lại (concatenate) thành một vector dài hơn — nhưng cách này có nhược điểm: 3 phần chỉ "nằm cạnh nhau", không có tương tác qua lại giữa chúng.

**Giải pháp — Kronecker product (tích Kronecker)**: đây là một phép nhân giữa 2 vector (hoặc nhiều vector) sao cho **kết quả chứa tất cả các cặp tích của từng phần tử** — nói cách khác, nó "nhân chéo" mọi thành phần của vector này với mọi thành phần của vector kia, tạo ra một không gian lớn hơn nhiều nhưng nắm bắt được **mọi tương tác có thể có** giữa các nguồn thông tin.

**Ví dụ minh họa Kronecker product với vector rất nhỏ (để hiểu cơ chế, không phải số thật của mô hình)**: giả sử a = [a1, a2] (2 chiều) và b = [b1, b2, b3] (3 chiều). Tích Kronecker a ⊗ b là:
```
a ⊗ b = [a1·b1, a1·b2, a1·b3, a2·b1, a2·b2, a2·b3]
```
Kết quả có 2×3=6 chiều — mỗi phần tử là tích của một phần tử từ a với một phần tử từ b, **bao trọn mọi cặp tương tác có thể có**. Với số cụ thể a=[2, 3], b=[1, 0, 4]:
```
a ⊗ b = [2×1, 2×0, 2×4, 3×1, 3×0, 3×4] = [2, 0, 8, 3, 0, 12]
```

**Vấn đề nảy sinh**: khi Kronecker 3 vector với nhau (x'_v, o_τv,i, ω_v), số chiều kết quả bằng **tích số chiều của cả 3 vector** — ví dụ nếu mỗi vector 100 chiều thì kết quả có 100×100×100 = 1 triệu chiều. Đây là bùng nổ tổ hợp (combinatorial explosion), không thể dùng trực tiếp.

**Giải pháp cho vấn đề này — Low-rank decomposition (nén hạng thấp, lấy cảm hứng từ LoRA)**: thay vì tính tường minh toàn bộ vector 1 triệu chiều rồi mới nén, ta **nén ngay từ đầu** bằng cách học các "phiên bản rút gọn" của từng vector con trước, rồi mới Kronecker các phiên bản rút gọn đó — tránh phải tính con số khổng lồ giữa chừng.

**Công thức**:
```
H_v = x'_v ⊗ o_τv,i ⊗ ω_v                     (Kronecker 3 vector — về mặt lý thuyết)
h_v = Σ_{i=1}^{r} ( w_i^(x)·x'_v ) ⊗ ( w_i^(o)·o_τv,i ) ⊗ ( w_i^(ω)·ω_v )   (tính thực tế, tránh bùng nổ)
```
`r` là "hạng" (rank) của phép nén — một số nhỏ (thực nghiệm gốc của BG-HGNN cho thấy r=4 hoặc r=5 là đủ). `w_i^(x)`, `w_i^(o)`, `w_i^(ω)` là các trọng số **có học được** (đây là những tham số thật sự cần train trong toàn bộ Nhánh A, số lượng tham số này rất nhỏ và không phụ thuộc vào số loại quan hệ).

**Trực giác đơn giản của công thức thứ hai**: thay vì làm phép Kronecker một lần với vector đầy đủ (rất tốn), ta làm r lần phép Kronecker với các **phiên bản đã co gọn** (rẻ hơn nhiều lần), rồi cộng r kết quả đó lại. Đây là kỹ thuật toán học chuẩn để xấp xỉ một phép tính "đắt" bằng tổng của nhiều phép tính "rẻ hơn nhiều" mà không mất quá nhiều độ chính xác — ý tưởng giống hệt LoRA (Low-Rank Adaptation) đang dùng phổ biến để fine-tune các mô hình ngôn ngữ lớn tiết kiệm tham số.

**Kết quả cuối bước này**: mỗi node v có một vector h_v duy nhất, đã chứa đựng thông tin từ (feature gốc, loại thực thể, bối cảnh quan hệ có trọng số bằng chứng) hòa trộn với nhau.

### 2.5. Lan truyền qua GNN nền — điểm chèn quan trọng thứ hai

**Vấn đề cần giải quyết**: h_v ở bước trên mới chỉ là đặc trưng "riêng" của từng node, tính từ hàng xóm trực tiếp một cách gộp chung (qua ω_v). Nhưng một GNN thật sự cần **lan truyền nhiều bước (layer)**, để thông tin từ hàng xóm-của-hàng-xóm cũng dần dần "lan tới" node v. Bước này còn là nơi thứ hai có thể chèn trọng số bằng chứng, và **quan trọng không kém bước 2.3** — vì bước 2.3 chỉ tính vector "bối cảnh" một lần trước khi vào GNN, còn bước 2.5 quyết định cách thông tin thực sự lan truyền qua nhiều bước trong suốt quá trình học.

**Vì sao dùng GAT (Graph Attention Network) làm nền**: GAT là một loại GNN mà khi một node tổng hợp thông tin từ hàng xóm, nó **tự học ra một trọng số attention** cho từng hàng xóm (hàng xóm nào "quan trọng hơn" thì được nghe nhiều hơn) — cơ chế này rất tự nhiên để chèn thêm trọng số bằng chứng vào cùng lúc.

**Cách tính attention gốc của GAT (không có trọng số bằng chứng)**:
```
e_uv = LeakyReLU( a^T [W h_u ‖ W h_v] )
α_uv = softmax_u(e_uv)
h_v^(mới) = σ( Σ_u α_uv · W h_u )
```
Giải thích từng dòng:
- Dòng 1: với mỗi cặp (u,v), tính một "điểm thô" e_uv đo mức độ hàng xóm u nên được chú ý bởi v. `W` là ma trận trọng số học được (biến đổi tuyến tính), `a` là một vector trọng số học được khác, `‖` là phép nối 2 vector lại, LeakyReLU là hàm phi tuyến thông thường.
- Dòng 2: **softmax** biến các điểm thô e_uv (của mọi hàng xóm u của v) thành một phân phối xác suất — tức là các số dương, cộng lại bằng 1. Hàng xóm nào có điểm thô cao hơn thì sau softmax có tỷ trọng α_uv lớn hơn.
- Dòng 3: node v cập nhật giá trị mới bằng cách cộng có trọng số (theo α_uv) các giá trị W h_u của mọi hàng xóm.

**Vấn đề**: cơ chế attention này tự học hoàn toàn từ dữ liệu, không biết gì về "cạnh này có bằng chứng lâm sàng, cạnh kia chỉ là suy luận" — hai cạnh có thể được gán trọng số attention y hệt nhau nếu đặc trưng của chúng "trông giống nhau" với mô hình.

**Cách chèn trọng số bằng chứng**:
```
e'_uv = e_uv + log(w_uv + ε)
α_uv = softmax_u(e'_uv)
```
Chỉ thêm một số hạng `log(w_uv + ε)` vào điểm thô trước khi đưa qua softmax.

**Vì sao cộng log(w) lại có tác dụng như nhân trọng số** (đây là điểm hay ho đáng giải thích kỹ): tính chất toán học của hàm mũ là `exp(e + log w) = exp(e) × exp(log w) = exp(e) × w`. Vì softmax về bản chất là chuẩn hóa các giá trị exp(e_uv), nên khi ta cộng log(w_uv) vào e_uv trước khi đưa vào softmax, kết quả sau softmax **tỷ lệ thuận trực tiếp** với w_uv — nói cách khác, hàng xóm có w_uv = 2 sẽ có tỷ trọng gấp đôi hàng xóm có w_uv = 1 (nếu điểm thô e_uv ban đầu bằng nhau), đúng như ta mong muốn, mà vẫn giữ được tính chất "tổng các tỷ trọng bằng 1" của softmax.

**Ví dụ số cụ thể**: giả sử node v có 2 hàng xóm, điểm thô GAT tính ra bằng nhau (e_u1v = e_u2v = 1.0, tức là GAT "không thiên vị" ai cả nếu chỉ nhìn đặc trưng), nhưng w_u1v = 1.0 (lâm sàng), w_u2v = 0.25 (suy luận).

Không chèn trọng số:
```
α_u1v = exp(1.0) / (exp(1.0) + exp(1.0)) = 0.5
α_u2v = exp(1.0) / (exp(1.0) + exp(1.0)) = 0.5
```
Hai hàng xóm được nghe như nhau — dù một cái đáng tin hơn hẳn.

Có chèn trọng số (ε bỏ qua vì không ảnh hưởng đáng kể ở đây):
```
e'_u1v = 1.0 + log(1.0)  = 1.0 + 0      = 1.0
e'_u2v = 1.0 + log(0.25) = 1.0 - 1.386  = -0.386

α_u1v = exp(1.0) / (exp(1.0)+exp(-0.386)) = 2.718 / (2.718+0.680) = 2.718/3.398 = 0.800
α_u2v = exp(-0.386) / (exp(1.0)+exp(-0.386)) = 0.680/3.398 = 0.200
```
Sau khi chèn, hàng xóm có bằng chứng mạnh (u1) được nghe với tỷ trọng 0.80 thay vì 0.5, hàng xóm bằng chứng yếu (u2) chỉ còn 0.20 — đúng theo tỷ lệ 1.0 : 0.25 = 4:1 (0.80/0.20 = 4). Đây là bằng chứng cho thấy cơ chế `log(w)` hoạt động chính xác như một phép nhân trọng số ẩn bên trong softmax.

**Chạy L lớp (layer)**: lặp lại toàn bộ bước 2.5 này L lần (khuyến nghị thử L=2 hoặc L=3 trước, theo dõi hiện tượng "over-smoothing" — càng nhiều lớp, các node càng có xu hướng trở nên giống nhau, làm mất đi sự phân biệt, hiện tượng này đã được HGHDA và HTINet2 ghi nhận khi tăng số lớp).

**Kết quả cuối Nhánh A**: mỗi node (herb, disease, và cả các loại thực thể khác) có một embedding cuối cùng, ký hiệu E_h^BG (cho herb) và E_d^BG (cho disease), đã mang thông tin lan truyền từ toàn bộ đồ thị, có ưu tiên theo bằng chứng.

> **Ablation bắt buộc để kiểm chứng**: chạy 3 biến thể — (a) chỉ trọng số hóa ω_v (mục 2.3), giữ GAT không trọng số (mục 2.5 dùng công thức gốc); (b) chỉ trọng số hóa GAT, giữ ω_v trung bình đều; (c) cả hai cùng lúc (bản đầy đủ) — để biết cơ chế nào đóng góp thật sự, tránh trường hợp chèn trọng số ở 2 nơi nhưng thực ra chỉ 1 nơi có tác dụng (hoặc tệ hơn, "double-count" khiến hiệu ứng bị phóng đại quá mức).

---

## 3. NHÁNH B — Cross-channel Attention (diễn giải phân tử)

Nhánh này chạy **song song** với Nhánh A, nhưng chỉ tập trung vào đúng một cặp (herb, disease) đang được xét tại một thời điểm — không xử lý toàn bộ đồ thị.

### 3.1. Trích subset liên quan

**Ý tưởng**: với một herb h cụ thể (ví dụ "Đương quy"), ta lấy ra đúng các ingredient mà nó chứa (ví dụ Ferulic acid, Ligustilide...). Với một disease d cụ thể (ví dụ "Alzheimer"), ta lấy ra đúng các target protein liên quan tới nó (ví dụ APOE, TLR4...). Đây là bước "khoanh vùng" — chỉ nhìn vào đúng những thực thể có liên hệ trực tiếp, thay vì toàn bộ đồ thị.

**Nguồn dữ liệu cho bước này**: thay vì phải tính lại embedding riêng cho ingredient/target (như iCAM-Net gốc phải làm bằng autoencoder + hypergraph riêng), ta **tái sử dụng luôn embedding đã có từ Nhánh A** (vì ingredient và target vốn dĩ cũng là node trong KG 9 thực thể, nên đã có embedding từ bước 2.5). Đây là điểm tiết kiệm tính toán so với thiết kế gốc của iCAM-Net.

Ký hiệu: E_C^h ∈ R^(N_c × d) là ma trận embedding của N_c ingredient thuộc herb h (mỗi hàng là embedding d chiều của một ingredient, lấy từ Nhánh A). E_P^d ∈ R^(N_p × d) tương tự cho N_p target thuộc disease d.

### 3.2. Attention hai chiều Q/K/V (Query/Key/Value)

**Trực giác về cơ chế Q/K/V** (giải thích tổng quát trước khi vào công thức): đây là cơ chế giống hệt attention trong kiến trúc Transformer. Ý tưởng là: mỗi ingredient "hỏi" (Query) tất cả các target "chìa khóa của bạn (Key) khớp với câu hỏi của tôi tới mức nào?", mức khớp này quyết định ingredient đó nên "nghe" giá trị (Value) của target nào nhiều hơn. Quan hệ này chạy theo cả 2 chiều: ingredient hỏi target, và ngược lại target cũng hỏi ingredient.

**Công thức**:
```
Q_C, K_C, V_C = Linear(E_C^h)        Q_P, K_P, V_P = Linear(E_P^d)

F_C = softmax( Q_C K_P^T / √d_attn ) V_P     (chiều: ingredient hỏi target)
F_P = softmax( Q_P K_C^T / √d_attn ) V_C     (chiều: target hỏi ingredient)
```
- `Linear(...)` là các phép biến đổi tuyến tính riêng biệt (có tham số học được, nhưng số lượng tham số này cố định, không phụ thuộc vào bao nhiêu ingredient/target).
- `Q_C K_P^T` là tích ma trận, cho ra một ma trận điểm số kích thước N_c × N_p — mỗi phần tử (i,j) là "điểm khớp" giữa ingredient thứ i và target thứ j.
- Chia cho `√d_attn` là một bước chuẩn hóa thường thấy trong Transformer, giúp giá trị không quá lớn trước khi vào softmax (tránh gradient bị bão hòa).
- `softmax(...)` áp dụng theo hàng (mỗi ingredient), biến điểm số thành tỷ trọng.
- Nhân với `V_P` để tổng hợp thông tin từ target theo tỷ trọng đó, ra F_C — embedding của ingredient đã "được làm giàu" bằng thông tin từ các target liên quan.

**Kết quả phụ quan trọng**: ma trận A_cp (attention ingredient→target) và A_pc (attention target→ingredient) — đây chính là "bằng chứng số" cho việc ingredient nào tương tác mạnh với target nào, dùng để diễn giải sau này (giống case study 5′-GMP–APOE của iCAM-Net).

### 3.3. Attention-driven pooling — gộp thành 1 vector cho herb/disease

**Vấn đề**: sau bước 3.2, ta có F_C (N_c vector, một cho mỗi ingredient) chứ chưa có **một** vector duy nhất đại diện cho cả herb h. Cần "gộp" N_c vector này lại.

**Cách gộp**: không gộp đều (trung bình đơn giản), mà ưu tiên ingredient nào có tương tác mạnh nhất với ít nhất một target nào đó (lấy giá trị lớn nhất theo hàng của ma trận attention, rồi chuẩn hóa lại bằng softmax).

```
a_C = softmax( max(A_cp, dim=1) )
E_h^attn = Σ_k (a_C)_k · (F_C)_k
```
Đọc thành lời: với mỗi ingredient k, tìm giá trị attention lớn nhất mà nó có với bất kỳ target nào (max theo dòng) — đây là "độ nổi bật" của ingredient đó. Chuẩn hóa các độ nổi bật này thành tỷ trọng bằng softmax, rồi dùng tỷ trọng đó để cộng có trọng số các vector F_C, ra một vector E_h^attn duy nhất đại diện cho herb h dưới góc nhìn "tương tác phân tử với disease d cụ thể này".

Làm tương tự để có E_d^attn cho disease d (dùng A_pc và F_P).

**Kết quả cuối Nhánh B**: E_h^attn, E_d^attn (mỗi cặp herb-disease cụ thể sẽ cho ra một cặp embedding riêng, vì nhánh này chạy theo từng cặp chứ không phải một lần cho toàn đồ thị như Nhánh A).

---

## 4. NHÁNH C — Kết hợp hai nhánh (Fusion)

**Vấn đề**: giờ có 2 loại thông tin về cùng một herb h — E_h^BG (từ Nhánh A, mang bối cảnh toàn đồ thị: symptom, meridian, prescription...) và E_h^attn (từ Nhánh B, mang tương tác phân tử cụ thể với disease d đang xét). Cần gộp lại thành một vector cuối cùng.

**Cách 1 — Concat + MLP (mạng nhiều lớp nhỏ)**:
```
E_h^final = MLP_fuse( [ E_h^BG ‖ E_h^attn ] )
```
Nối 2 vector lại, rồi cho qua một mạng neural nhỏ để "học cách trộn" 2 nguồn thông tin thành một biểu diễn thống nhất.

**Cách 2 — Gated sum (cổng trọng số học được)**, có thể thử nếu Cách 1 không ổn định:
```
g = sigmoid( W_g [E_h^BG ‖ E_h^attn] )
E_h^final = g ⊙ E_h^BG + (1-g) ⊙ E_h^attn
```
Ở đây `g` là một vector "cổng" (giá trị từ 0 đến 1 theo từng chiều, nhờ hàm sigmoid), quyết định ở mỗi chiều nên tin theo E_h^BG hay E_h^attn nhiều hơn — mô hình tự học ra cổng này. `⊙` là phép nhân từng phần tử (element-wise).

Làm tương tự cho disease để ra E_d^final.

> **Câu hỏi thiết kế cần thử nghiệm**: nên lấy E_h^BG từ layer nào của Nhánh A để đưa vào đây? Có 2 lựa chọn: (a) lấy ngay sau bước 2.4 (Kronecker + nén hạng thấp, "gần gốc", còn giữ nguyên tín hiệu phân tử thô, chưa bị pha loãng qua nhiều bước lan truyền), hoặc (b) lấy sau toàn bộ L lớp GAT ở bước 2.5 ("xa gốc", đã có bối cảnh toàn đồ thị phong phú hơn nhưng có thể bị "làm mờ" thông tin chi tiết). Khuyến nghị thử cả hai và so sánh bằng thực nghiệm.

---

## 5. NHÁNH D — Multi-task heads và hàm mất mát (loss)

### 5.1. Task chính — dự đoán Herb-Disease Association

**Công thức**:
```
F_HD = sigmoid( E_h^final · (E_d^final)^T )
```
`E_h^final · (E_d^final)^T` là tích vô hướng (dot product) giữa 2 vector — một con số duy nhất, càng lớn nghĩa là 2 vector càng "giống hướng nhau" (tương đồng cao). `sigmoid(...)` ép con số này về khoảng (0,1), để đọc như một xác suất.

**Hàm mất mát (Binary Cross-Entropy — BCE)**: đây là hàm đo "mức độ sai" giữa xác suất dự đoán F_HD và nhãn thật y_HD (y_HD=1 nếu thực sự có association, =0 nếu không). Trực giác: nếu mô hình dự đoán xác suất gần đúng với nhãn thật, loss sẽ nhỏ; nếu dự đoán sai lệch nhiều (ví dụ dự đoán 0.1 nhưng thực ra có association), loss sẽ rất lớn — hàm này phạt nặng những dự đoán "tự tin nhưng sai".
```
L_HD = -[ y_HD·log(F_HD) + (1-y_HD)·log(1-F_HD) ]     (trung bình trên cả batch)
```

### 5.2. Task phụ — dự đoán Ingredient-Target Association (đóng vai trò như CPA của iCAM-Net)

**Vì sao cần task phụ này**: nếu chỉ train theo task chính (HDA), mô hình có thể đạt điểm số dự đoán tốt mà **không cần** học các mối tương tác ingredient-target có ý nghĩa thật sự — nó có thể "gian lận" bằng cách tìm ra các đường tắt thống kê nào đó trong dữ liệu. Việc thêm một task phụ, buộc mô hình phải dự đoán đúng luôn cả association ở mức phân tử (ingredient-target), ép nhánh attention (Nhánh B) phải học ra các trọng số có ý nghĩa sinh học thật, chứ không chỉ là con số ngẫu nhiên giúp task chính đạt điểm cao.

**Công thức**: dùng trực tiếp F_C, F_P đã tính ở bước 3.2:
```
F_IT = sigmoid( Bilinear(F_C, F_P) )
L_IT = BCE(F_IT, y_IT)
```
`Bilinear(a, b) = a^T · M · b` (M là ma trận tham số học được) — một cách tổng quát hơn dot product thường, cho phép mô hình học ra cách "a và b tương tác" phức tạp hơn là chỉ nhân đơn giản.

### 5.3. Hàm mất mát tổng và cách huấn luyện xen kẽ

```
L_total = L_HD + α · L_IT
```
`α` là một số cân bằng giữa 2 task (task phụ không nên "lấn át" task chính) — cần dò tìm bằng thực nghiệm, thử α ∈ {0.01, 0.05, 0.1, 0.2, 0.3}, không mặc định lấy 0.05 của iCAM-Net vì cấu trúc KG khác hẳn.

**Huấn luyện xen kẽ (alternating optimization)**: vì số lượng mẫu của 2 task khác nhau rất nhiều, cách train ổn định là: mỗi epoch, chạy hết một vòng batch của task HDA (cập nhật tham số), rồi chạy hết một vòng batch của task IT (cập nhật tiếp tham số), lặp lại — thay vì trộn lẫn 2 loại batch cùng lúc.

---

## 6. Kế hoạch huấn luyện, đánh giá và ablation chi tiết

### 6.1. Negative sampling
Chọn ngẫu nhiên các cặp herb-disease không có trong nhãn dương, số lượng bằng đúng số mẫu dương — đây là cách 3/4 bài đã khảo sát (HTINet2, HGHDA, iCAM-Net) đồng thuận sử dụng.

### 6.2. Ba kiểu Cross-validation
- **Global 5-fold**: chia ngẫu nhiên toàn bộ cặp herb-disease thành 5 phần, mỗi lần lấy 1 phần làm test, 4 phần còn lại train — đo khả năng dự đoán khi cả herb và disease đều đã "quen mặt" lúc train.
- **Local 5-fold theo herb**: chia herb thành 5 nhóm, mỗi lần chọn 1 nhóm herb, giữ lại toàn bộ association của các herb đó làm test, train trên phần còn lại — đo khả năng dự đoán khi có herb hoàn toàn mới, chưa từng thấy lúc train.
- **Local 5-fold theo disease**: tương tự nhưng chia theo disease.
- **Cold-start test**: bỏ hẳn 1 herb cụ thể (toàn bộ association của nó) ra khỏi tập train, sau đó xem mô hình dự đoán lại các association đó tốt tới đâu — so sánh với tỷ lệ khôi phục ~88-90% mà HGHDA và iCAM-Net đã đạt được.

### 6.3. Kiểm tra tính robust bằng xáo trộn ngẫu nhiên (theo mẫu HGHDA — đa mức, mạnh hơn kiểu xáo trộn 1 lần của iCAM-Net)
Chọn ngẫu nhiên và tráo đổi 10%, 20%, ..., 90% các phần tử trong ma trận kề herb-disease, rồi chạy lại 5-fold CV ở từng mức xáo trộn. Kỳ vọng: AUROC/AUPRC giảm dần và tiệm cận 0.5 (mức ngẫu nhiên hoàn toàn) khi xáo trộn tới 90% — nếu đúng như vậy, chứng minh được mô hình thực sự học từ cấu trúc sinh học thật trong dữ liệu, chứ không phải chỉ ghi nhớ (overfit) các mẫu ngẫu nhiên.

### 6.4. Bảng ablation đầy đủ (mỗi dòng chỉ đổi một thứ so với bản đầy đủ, để tách bạch từng đóng góp)

| Biến thể | Thay đổi so với bản đầy đủ | Mục đích đo lường |
|---|---|---|
| Full model | (không đổi) | baseline để so sánh |
| − trọng số ở ω_v | dùng công thức trung bình đều (mục 2.3, công thức gốc) | đóng góp của trọng số bằng chứng ở bước gộp bối cảnh quan hệ |
| − trọng số ở GAT | bỏ số hạng log(w_uv) (mục 2.5, công thức gốc) | đóng góp của trọng số bằng chứng ở bước lan truyền |
| − nhánh B hoàn toàn | chỉ dùng E^BG, bỏ qua E^attn ở bước Fusion | đóng góp của nhánh attention/khả năng diễn giải |
| − nhánh A, thay bằng feature khởi tạo ngẫu nhiên | random init thay vì BG-HGNN | đóng góp của bối cảnh toàn KG (so sánh với thí nghiệm HGHDA-s) |
| − multi-task (đặt α=0) | chỉ tối ưu L_HD | vai trò của supervision phân tử (so sánh với iCAM-m) |
| xáo trộn thật mapping herb-ingredient/disease-target | tráo ngẫu nhiên nhưng giữ nguyên số lượng | cấu trúc sinh học thật có thực sự cần thiết không (so sánh với iCAM-r) |
| one-hot thay dense random projection | đổi cách mã hóa loại thực thể/quan hệ ở mục 2.2 | xác nhận lại kết luận gốc của BG-HGNN có đúng trên KG mới hay không |

### 6.5. Baseline so sánh
HGHDA, iCAM-Net (cài đặt lại/adapt trên HERB 2.0), HDAPM-NCP, và R-GCN thuần (không có trọng số bằng chứng) — baseline R-GCN thuần đặc biệt quan trọng để làm rõ: cải thiện đến từ việc **có trọng số bằng chứng** hay đến từ **việc đổi sang backbone BG-HGNN**, hai điều này cần tách biệt rõ trong phần kết quả.

---

## 7. Rủi ro và điểm cần theo dõi khi triển khai thực tế

1. **Chi phí tính toán cộng dồn**: Nhánh B vẫn có độ phức tạp O(m×n) theo từng cặp herb-disease (đây là giới hạn mà chính iCAM-Net tự thừa nhận trong phần Limitations của họ) — cần đo thời gian chạy thực tế của toàn bộ pipeline (Nhánh A + B + C + D cộng lại), không chỉ dựa vào con số throughput mà BG-HGNN tự báo cáo (họ chỉ đo riêng backbone, chưa cộng thêm nhánh attention).
2. **Trùng lặp thông tin giữa 2 nhánh**: vì Nhánh A cũng đã lan truyền qua quan hệ herb-ingredient (nó là 1 trong 28 quan hệ), E_h^BG có thể đã gián tiếp "biết" một phần thông tin ingredient rồi. Ablation "bỏ nhánh B hoàn toàn" (mục 6.4) là bắt buộc để chứng minh nhánh B đóng góp thêm thông tin thật sự, không chỉ lặp lại cái Nhánh A đã có.
3. **Cách chọn giá trị w_uv cho 4 tầng bằng chứng**: cần một bảng ánh xạ rõ ràng và có lý do giải thích (không chỉ chọn số tùy ý) — có thể tham khảo cách HTINet2 tính "inferred evidence score" (dùng p-value + ngưỡng trung bình) làm cơ sở tham khảo, dù cách tiếp cận của họ là số liên tục từ một công thức, còn ở đây là ánh xạ rời rạc theo 4 tầng.
4. **Nguy cơ double-count trọng số**: vì có 2 điểm chèn trọng số (mục 2.3 và 2.5), cần ablation rõ ràng (đã có trong bảng mục 6.4) trước khi kết luận trọng số bằng chứng thực sự giúp ích, tránh trường hợp 2 nơi cùng khuếch đại một hiệu ứng khiến mô hình bị lệch quá mức về phía các cạnh có bằng chứng mạnh, bỏ qua hoàn toàn tín hiệu từ các cạnh bằng chứng yếu.
5. **Siêu tham số r (rank nén) và L (số lớp GAT)**: cần dò tìm riêng trên KG HERB 2.0, không mặc định theo giá trị gốc r=4-5, L=3 của BG-HGNN/HGHDA — vì đây là dataset khác hẳn, phân phối số lượng quan hệ trên mỗi node cũng khác, giá trị tối ưu có thể khác.
